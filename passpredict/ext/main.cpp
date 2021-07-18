#include <iostream>
#include <math.h>

#include "sofa.h"

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
        // constructor
        Location(double alat, double alon, double ah){
            lat = alat;
            lon = alon;
            h = ah;
            return;
        }
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




int main(){
    int i;
    std::cout << "hello\n";

    // Define TLE as two strings

    // Put TLE orbit data into appropriate cpp data structure

    // Put Location information into appropriate cpp data structure

    // Put TLE and Location data structures into Observer cpp data structure

    // Propagate satellite

    // Find location ecef position

    // Find az, el, range of observer

    return 0;
}