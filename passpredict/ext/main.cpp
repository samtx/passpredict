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
    int j;
    double u1, u2, tt1, tt2;
    j = iauDtf2d ( "UTC", 2010, 7, 24, 11, 18, 7.318, &u1, &u2 );
    j = passpredict::Utc2tt(u1, u2, tt1, tt2);

    std::cout << std::fixed << std::setprecision(4) << "u1 = " << u1 << "   u2 = " << std::setprecision(12)  << u2 << std::endl;
    std::cout << std::fixed << std::setprecision(4) << "tt1 = " << tt1 << "   tt2 = " << std::setprecision(12)  << tt2 << std::endl;


    lat = 30.1990;
    lon = -97.8616;
    h = 0.0;
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

    // // Put TLE and Location data structures into Observer cpp data structure
    // // ISS (ZARYA)
    // char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    // char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    // passpredict::Orbit sat (tle1, tle2);

    // const double deg2rad = PASSPREDICT_DEG2RAD;
	// const double xpdotp = 1440.0 / PASSPREDICT_2PI;

    // Find AOS


    char tle1[] = "1 25544U 98067A   21215.85360061  .00001804  00000-0  41053-4 0  9990";
    char tle2[] = "2 25544  51.6430 107.0086 0001317 252.4656 228.3182 15.48870648295985";
    passpredict::Orbit orbit (tle1, tle2);
    passpredict::Satellite satellite (orbit);
    passpredict::Observer observer(location, satellite);
    int year = 2021, mon = 8, day = 3, hr = 0, minute = 0;
    double sec = 0.0;
    double jd, jdFrac, aos;
    double tmax;
    SGP4Funcs::jday_SGP4(year, mon, day, hr, minute, sec, jd, jdFrac);
    tmax = jd + 10;
    std::cout << "jd=" << jd << "  tmax=" << tmax << std::endl;
    aos = passpredict::FindAOS(jd, tmax, location, satellite);
    std::cout << "AOS = " << aos << std::endl;
    std::cout << "... Predict ... " << std::endl;

    // std::list<std::shared_ptr<passpredict::Overpass>> overpasses;
    // std::shared_ptr<passpredict::Overpass> overpass;
    // overpasses = passpredict::Predict(location, satellite, jd, tmax);
    // int n_overpasses = overpasses.size();

    // std::cout << "Found " << n_overpasses << " overpasses" << std::endl;
    // for (auto const& overpass : overpasses) {
    //     overpass->PrintLn();
    // }


    return 0;
}