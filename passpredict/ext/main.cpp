#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>

extern "C"
{
#include "sofa.h"
};

#include "SGP4.h"
#include "passpredict.h"


#define D_R_EARTH 6378.137        //  Earth mean equatorial radius [km]
#define D_e2_EARTH 0.006694385000 // Earth eccentricity squared

struct Omm {
    char satnum[9];
    double jdsatepoch;
    double jdsatepochF;
    double bstar;
    double inclo;  // deg
    double nodeo;  // deg
    double ecco;
    double argpo;  // deg
    double mo;     // deg
    double no_kozai;
    int revnum;
    int elnum;
    char classification;
    int ephtype;
};

class Orbit
{
private:
    gravconsttype m_whichconst = wgs84;

public:
    const char *tle1;
    const char *tle2;
    elsetrec satrec;

    Orbit() {};

    Orbit(const Omm &omm)
    {
        bool err;
        double ndot, nddot;  // These aren't used in sgp4 so they are dummy variables

        const double deg2rad = PASSPREDICT_DEG2RAD;
        const double xpdotp = 1440.0 / PASSPREDICT_2PI;

        strcpy(satrec.satnum, omm.satnum);
        satrec.jdsatepoch = omm.jdsatepoch;
        satrec.jdsatepochF = omm.jdsatepochF;
		satrec.bstar = omm.bstar;
        satrec.ecco = omm.ecco;
        satrec.argpo = omm.argpo * deg2rad;
        satrec.inclo = omm.inclo * deg2rad,
        satrec.mo = omm.mo * deg2rad;
        satrec.no_kozai = omm.no_kozai / xpdotp;
        satrec.nodeo = omm.nodeo * deg2rad;
        satrec.revnum = omm.revnum;

        err = SGP4Funcs::sgp4init(
            m_whichconst, 'i', satrec.satnum, satrec.jdsatepoch + satrec.jdsatepochF - 2433281.5,
		    satrec.bstar, ndot, nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo,
            satrec.no_kozai, satrec.nodeo, satrec
		);
    };

    Orbit(char* atle1, char* atle2)
    {
        double dummy;
        SGP4Funcs::twoline2rv(
            atle1, atle2, ' ', ' ', 'i', m_whichconst, dummy, dummy, dummy, satrec
        );
    };
};

class Satellite
{
private:
    double m_rteme[3] = {0, 0, 0};  // TEME position vector at time jd
    double m_vteme[3] = {0, 0, 0};  // TEME velocity vector at time jd
    double m_jd;        // julian date
    double m_tsince;    // minutes since epoch
    double m_epoch;     // epoch, julian date
public:
    Orbit orbit;    // orbit data, OMM/TLE
    std::string name; // satellite name
    int satid;        // NORAD satellite ID number
    double recef[3];  // ECEF position vector at time jd
    double vecef[3];  // ECEF velocity vector at time jd

    Satellite(Orbit aorbit)
    {
        orbit = aorbit;
        m_epoch = orbit.satrec.jdsatepoch + orbit.satrec.jdsatepochF;
    };

    void propagate_tsince(double t_tsince)
    {
        // find time since julian date in minutes
        m_jd = m_epoch + (t_tsince / 1440.0);
        m_tsince = t_tsince;
        SGP4Funcs::sgp4(orbit.satrec, m_tsince, m_rteme, m_vteme);
    }

    /*
    Propagate the Satellite to time t

    inputs:
        t       time in julian date
    */
    void propagate_jd(double t_jd)
    {
        // find time since julian date in minutes
        m_jd = t_jd;
        m_tsince = (m_epoch - m_jd) / 1440.0;
        SGP4Funcs::sgp4(orbit.satrec, m_tsince, m_rteme, m_vteme);
    }


    void print()
    {
        int i;
        using namespace std;
        cout.precision(10);
        cout << "jd    =  " << fixed << m_jd << endl;
        cout.precision(4);
        cout << "rteme = [" << fixed;
        for (i=0; i<3; i++) {
            cout << m_rteme[i] << ", ";
        }
        cout << "]" << endl;
        cout << "vteme = [" << fixed;
        for (i=0; i<3; i++)
            cout << m_vteme[i] << ", ";
        cout << "]" << endl;
    };

    void print_oneline()
    {
        int i;
        using namespace std;
        cout.precision(8);
        cout << "   " << setw(14) << m_tsince;
        for (i=0; i<3; i++)
            cout << "   " << setw(14) << m_rteme[i];
        for (i=0; i<3; i++)
            cout << "   " << setw(12) << m_vteme[i];
    }
};

class Location
{
public:
    double lat;      // latitude (degrees)
    double lon;      // longitude (degrees)
    double h;        // height above MSL (m)
    double recef[3]; // ECEF position vector (km)

    // constructor
    Location(double alat, double alon, double ah)
    {
        lat = alat;
        lon = alon;
        h = ah;
        return;
    };

    // get ECEF position vector
    void site_ecef()
    {
        //  Vallado, Eq. 3-7
        double phi_gd_rad, h_ellp_km, C, S, r_delta, r_K, sin_phi_gd;
        double lmda_rad;
        phi_gd_rad = lat * DD2R;
        sin_phi_gd = std::sin(phi_gd_rad);
        h_ellp_km = h * 0.001; // convert height to km
        C = D_R_EARTH / std::sqrt(1 - D_e2_EARTH * std::pow(sin_phi_gd, 2));
        S = C * (1 - D_e2_EARTH);
        r_delta = (C + h_ellp_km) * std::cos(phi_gd_rad);
        r_K = (S + h_ellp_km) * sin_phi_gd;

        // Vallado, Alg 51, p.430
        lmda_rad = lon * DD2R;
        recef[0] = r_delta * std::cos(lmda_rad);
        recef[1] = r_delta * std::sin(lmda_rad);
        recef[2] = r_K;
    };
};

class Observer
{
public:
    double pos[3]; // position in ECI coordinates (km)
    double vel[3]; // velocity in ECI coordinates (km/s)
    double jd;     // julian date for time
    double el;     // elevation (degrees)
    double az;     // azimuth   (degrees)
    double range;  // range (km)
};

int main()
{
    int i;
    std::cout << "hello\n";

    // Location information
    double lat, lon, h;

    // Define TLE as two strings

    // Put TLE orbit data into appropriate cpp data structure

    // Put Location information into appropriate cpp data structure

    // phi_gd = 39.007  # [deg]
    //     lmda = -104.883  # [deg]
    //     alt = 2187.0  # [m]
    //     r_ECEF = topocentric.site_ECEF(phi_gd, lmda, alt)
    //     r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    //     for i in [0, 1, 2]:

    lat = 39.007;
    lon = -104.883;
    h = 2187.0;
    Location location(lat, lon, h);
    {
        using namespace std;
        cout << "Lat: " << location.lat << "\n";
        cout << "Lon: " << location.lon << "\n";
        cout << "H: " << location.h << "\n";
        // Find location ecef position
        location.site_ecef();
        cout << "recef: ";
        cout << fixed;
        for (i = 0; i < 3; i++)
        {
            cout << location.recef[i] << ", ";
        }
        cout << endl;
    }

    // Put TLE and Location data structures into Observer cpp data structure
    // ISS (ZARYA)
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    Orbit sat (tle1, tle2);

    const double pi = 3.14159265358979323846;
    const double deg2rad = PASSPREDICT_DEG2RAD;
	const double xpdotp = 1440.0 / PASSPREDICT_2PI;

    // Print satrec
    {
        using namespace std;
        cout << endl << "satrec from TLE strings" << endl;
        cout << "satrec.jdsatepoch = " << sat.satrec.jdsatepoch << endl;
        cout << "satrec.jdsatepochF = " << sat.satrec.jdsatepochF << endl;
        cout << "sgp4init epoch = " << (sat.satrec.jdsatepoch + sat.satrec.jdsatepochF) - 2433281.5 << endl;
        cout << "satrec.bstar = " << sat.satrec.bstar << endl;
        cout << "satrec.inclo = " << sat.satrec.inclo / deg2rad << endl;
        cout << "satrec.nodeo = " << sat.satrec.nodeo / deg2rad << endl;
        cout << "satrec.ecco = " << sat.satrec.ecco << endl;
        cout << "satrec.argpo = " << sat.satrec.argpo / deg2rad << endl;
        cout << "satrec.mo = " << sat.satrec.mo / deg2rad << endl;
        cout << "satrec.no_kozai = " << sat.satrec.no_kozai * xpdotp << endl;
        cout << "satrec.revnum = " << sat.satrec.revnum << endl;
    };

    // Use Omm
    Omm omm;
    strcpy(omm.satnum, "25544");        // satnum
    omm.jdsatepoch = 2.45942e+6;     // jdsatepoch
    omm.jdsatepochF = 0.469801;       // jdsatepochF
    omm.bstar = 4.2487e-5;      // bstar
    omm.inclo = 51.6426;        // inclo
    omm.nodeo = 178.1369;       // nodeo
    omm.ecco = 0.0001717;      // ecco
    omm.argpo = 174.7410;       // argpo
    omm.mo = 330.7918;       // mo
    omm.no_kozai = 15.4883;        // no_kozai
    omm.revnum = 293750;         // revnum
    omm.elnum = 993;            // elnum


    omm.classification = 'u';            // classification
    omm.ephtype = 0;               // ephtype

    Orbit sat2 (omm);
    // Print satrec
    {
        using namespace std;
        cout << endl << "satrec from Omm" << endl;
        cout << "satrec.jdsatepoch = " << sat2.satrec.jdsatepoch << endl;
        cout << "satrec.jdsatepochF = " << sat2.satrec.jdsatepochF << endl;
        cout << "sgp4init epoch = " << (sat2.satrec.jdsatepoch + sat2.satrec.jdsatepochF) - 2433281.5 << endl;
        cout << "satrec.bstar = " << sat2.satrec.bstar << endl;
        cout << "satrec.inclo = " << sat2.satrec.inclo / deg2rad << endl;
        cout << "satrec.nodeo = " << sat2.satrec.nodeo / deg2rad << endl;
        cout << "satrec.ecco = " << sat2.satrec.ecco << endl;
        cout << "satrec.argpo = " << sat2.satrec.argpo / deg2rad << endl;
        cout << "satrec.mo = " << sat2.satrec.mo / deg2rad << endl;
        cout << "satrec.no_kozai = " << sat2.satrec.no_kozai * xpdotp << endl;
        cout << "satrec.revnum = " << sat2.satrec.revnum << endl;
    };

    // Propagate satellite
    {
        char tle1[] = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
        char tle2[] = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
        double tsince;
        Orbit orbit (tle1, tle2);
        Satellite satellite (orbit);
        for (i=0; i<13; i++) {
            tsince = i * 360.0;
            satellite.propagate_tsince(tsince);
            satellite.print_oneline();
            std::cout << std::endl;
        }

    }



    // Find az, el, range of observer


    return 0;
}