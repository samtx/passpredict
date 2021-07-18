#include <iostream>
#include <math.h>

extern "C" {
    #include "sofa.h"
}

#include "sgp4.h"

#define D_R_EARTH 6378.137  //  Earth mean equatorial radius [km]
#define D_e2_EARTH 0.006694385000  // Earth eccentricity squared

typedef struct orbit_t {

} orbit_t;

class Satellite {
    public:
        orbit_t orbit;  // orbit data, OMM/TLE
        double epoch;   // epoch, julian date
        std::string name;    // satellite name
        int satid;  // NORAD satellite ID number
};


class Location {
    public:
        double lat;  // latitude (degrees)
        double lon;  // longitude (degrees)
        double h;    // height above MSL (m)
        double recef[3]; // ECEF position vector (km)
        // constructor
        Location(double alat, double alon, double ah){
            lat = alat;
            lon = alon;
            h = ah;
            return;
        };

        // get ECEF position vector
        void site_ecef() {
            //  Vallado, Eq. 3-7
            // phi_gd_rad = phi_gd * DEG2RAD  # convert [deg] to [rad]
            // h_ellp_km = h_ellp * 0.001  # convert [m] to [km]
            // C = R_EARTH / math.sqrt(1 - e2_EARTH * math.sin(phi_gd_rad) ** 2)
            // S = C * (1 - e2_EARTH)
            // r_delta = (C + h_ellp_km) * math.cos(phi_gd_rad)
            // r_K = (S + h_ellp_km) * math.sin(phi_gd_rad)
            // return (r_delta, r_K)
            double phi_gd_rad, h_ellp_km, C, S, r_delta, r_K, sin_phi_gd;
            double lmda_rad;
            phi_gd_rad = lat * DD2R;
            sin_phi_gd = std::sin(phi_gd_rad);
            h_ellp_km = h * 0.001;  // convert height to km
            C = D_R_EARTH / std::sqrt(1 - D_e2_EARTH * std::pow(sin_phi_gd,2));
            S = C * (1 - D_e2_EARTH);
            r_delta = (C + h_ellp_km) * std::cos(phi_gd_rad);
            r_K = (S + h_ellp_km) * sin_phi_gd;

            // Vallado, Alg 51, p.430
            /* 
            lmda_rad = math.radians(lmda)
            r_site_ECEF = np.array([r_delta * math.cos(lmda_rad), r_delta * math.sin(lmda_rad), r_K])
            return r_site_ECEF
            */
            lmda_rad = lon * DD2R;
            recef[0] = r_delta * std::cos(lmda_rad);
            recef[1] = r_delta * std::sin(lmda_rad);
            recef[2] = r_K;
        };
};




class Observer {
    public:
        double pos[3];  // position in ECI coordinates (km)
        double vel[3];  // velocity in ECI coordinates (km/s)
        double jd;  // julian date for time
        double el;  // elevation (degrees)
        double az;  // azimuth   (degrees)
        double range; // range (km)
}; 

/*
def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    phi_gd_rad = phi_gd * DEG2RAD  # convert [deg] to [rad]
    h_ellp_km = h_ellp * 0.001  # convert [m] to [km]
    C = R_EARTH / math.sqrt(1 - e2_EARTH * math.sin(phi_gd_rad) ** 2)
    S = C * (1 - e2_EARTH)
    r_delta = (C + h_ellp_km) * math.cos(phi_gd_rad)
    r_K = (S + h_ellp_km) * math.sin(phi_gd_rad)
    return (r_delta, r_K)

def site_ECEF(phi_gd, lmda, h_ellp):
    """Compute ECEF coordinates for tracking site on Earth

    References:
        Vallado, Algorithm 51, p.430
    """
    r_delta, r_K = site_declination_and_K(phi_gd, h_ellp)
    lmda_rad = math.radians(lmda)
    r_site_ECEF = np.array([r_delta * math.cos(lmda_rad), r_delta * math.sin(lmda_rad), r_K])
    return r_site_ECEF
*/


int main(){
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

    // Propagate satellite

    // Find location ecef position
    location.site_ecef();
    std::cout << "recef: ";
    std::cout << std::fixed;
    for (i=0; i<3; i++) {
        std::cout << location.recef[i] << ", ";
    }
    std::cout << std::endl;

    // Find az, el, range of observer

    return 0;
}