#include <iostream>
#include <math.h>

extern "C"
{
#include "sofa.h"
};

#include "SGP4.h"


#define D_R_EARTH 6378.137        //  Earth mean equatorial radius [km]
#define D_e2_EARTH 0.006694385000 // Earth eccentricity squared

class Orbit
{
private:
    elsetrec satrec;

public:
    const char *tle1;
    const char *tle2;

    Orbit(const char* atle1, const char* atle2)
    {
        double dummy;
        gravconsttype whichconst = wgs84;
        elsetrec asatrec;

        SGP4Funcs::twoline2rv(
            const_cast<char*>(atle1), const_cast<char*>(atle2),
            ' ', ' ', 'i', whichconst,
            dummy, dummy, dummy, asatrec
        );
    };
};

class Satellite
{
public:
    Orbit orbit;    // orbit data, OMM/TLE
    double epoch;     // epoch, julian date
    std::string name; // satellite name
    int satid;        // NORAD satellite ID number
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
    std::cout << "Lat: " << location.lat << "\n";
    std::cout << "Lon: " << location.lon << "\n";
    std::cout << "H: " << location.h << "\n";

    // Put TLE and Location data structures into Observer cpp data structure
    // ISS (ZARYA)
    const char* tle1 = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    const char* tle2 = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    Orbit sat (tle1, tle2);

    // Propagate satellite

    // Find location ecef position
    location.site_ecef();
    std::cout << "recef: ";
    std::cout << std::fixed;
    for (i = 0; i < 3; i++)
    {
        std::cout << location.recef[i] << ", ";
    }
    std::cout << std::endl;

    // Find az, el, range of observer


    return 0;
}