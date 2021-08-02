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
    passpredict::Location location(lat, lon, h);
    {
        using namespace std;
        cout << "Lat: " << location.lat_ << "\n";
        cout << "Lon: " << location.lon_ << "\n";
        cout << "H: " << location.h_ << "\n";
        // Find location ecef position
        cout << "recef: ";
        cout << fixed;
        for (i = 0; i < 3; i++)
        {
            cout << location.recef_[i] << ", ";
        }
        cout << endl;
    }

    // Put TLE and Location data structures into Observer cpp data structure
    // ISS (ZARYA)
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    passpredict::Orbit sat (tle1, tle2);

    const double deg2rad = PASSPREDICT_DEG2RAD;
	const double xpdotp = 1440.0 / PASSPREDICT_2PI;

    // Print satrec_
    {
        using namespace std;
        cout << endl << "satrec_ from TLE strings" << endl;
        cout << "satrec_.jdsatepoch = " << sat.satrec_.jdsatepoch << endl;
        cout << "satrec_.jdsatepochF = " << sat.satrec_.jdsatepochF << endl;
        cout << "sgp4init epoch = " << (sat.satrec_.jdsatepoch + sat.satrec_.jdsatepochF) - 2433281.5 << endl;
        cout << "satrec_.bstar = " << sat.satrec_.bstar << endl;
        cout << "satrec_.inclo = " << sat.satrec_.inclo / deg2rad << endl;
        cout << "satrec_.nodeo = " << sat.satrec_.nodeo / deg2rad << endl;
        cout << "satrec_.ecco = " << sat.satrec_.ecco << endl;
        cout << "satrec_.argpo = " << sat.satrec_.argpo / deg2rad << endl;
        cout << "satrec_.mo = " << sat.satrec_.mo / deg2rad << endl;
        cout << "satrec_.no_kozai = " << sat.satrec_.no_kozai * xpdotp << endl;
        cout << "satrec_.revnum = " << sat.satrec_.revnum << endl;
    };

    // Use Omm
    passpredict::Omm omm;
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

    passpredict::Orbit sat2 (omm);
    // Print satrec_
    {
        using namespace std;
        cout << endl << "satrec_ from Omm" << endl;
        cout << "satrec_.jdsatepoch = " << sat2.satrec_.jdsatepoch << endl;
        cout << "satrec_.jdsatepochF = " << sat2.satrec_.jdsatepochF << endl;
        cout << "sgp4init epoch = " << (sat2.satrec_.jdsatepoch + sat2.satrec_.jdsatepochF) - 2433281.5 << endl;
        cout << "satrec_.bstar = " << sat2.satrec_.bstar << endl;
        cout << "satrec_.inclo = " << sat2.satrec_.inclo / deg2rad << endl;
        cout << "satrec_.nodeo = " << sat2.satrec_.nodeo / deg2rad << endl;
        cout << "satrec_.ecco = " << sat2.satrec_.ecco << endl;
        cout << "satrec_.argpo = " << sat2.satrec_.argpo / deg2rad << endl;
        cout << "satrec_.mo = " << sat2.satrec_.mo / deg2rad << endl;
        cout << "satrec_.no_kozai = " << sat2.satrec_.no_kozai * xpdotp << endl;
        cout << "satrec_.revnum = " << sat2.satrec_.revnum << endl;
    };

    // Propagate satellite
    {
        char tle1[] = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
        char tle2[] = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
        double tsince;
        passpredict::Orbit orbit (tle1, tle2);
        passpredict::Satellite satellite (orbit);
        for (i=0; i<13; i++) {
            tsince = i * 360.0;
            satellite.PropagateTSince(tsince);
            satellite.PrintOneline();
            std::cout << std::endl;
        }

    }



    // Find az, el, range of observer


    return 0;
}