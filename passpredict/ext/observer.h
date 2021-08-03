#ifndef OBSERVER_H
#define OBSERVER_H

#include "location.h"
#include "satellite.h"
#include "passpredict.h"

namespace passpredict
{
    class Observer
    {
    public:
        Location location_;
        Satellite satellite_;
        double r_[3] = {0, 0, 0}; // position in ECI coordinates (km)
        double v_[3] = {0, 0, 0}; // velocity in ECI coordinates (km/s)
        double jd_;               // julian date for time
        double el_;               // elevation (degrees)
        double az_;               // azimuth   (degrees)
        double range_;            // range (km)

        Observer(Location location, Satellite satellite);
        double ComputeElevationAngle();
        void Ecef2Sez(double rho[3], double rsez[3]);
        void Sez2Razel(double rsez[3]);
        double FindLos(double, double);
        double FindAos(double, double);
        // void FindAos(double t0, double tmax);
        // void FindLos(double t0, double tmax);

    };
}
#endif