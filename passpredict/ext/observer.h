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
        Location location;
        Satellite satellite;
        double r[3] = {0, 0, 0}; // position in ECI coordinates (km)
        double v[3] = {0, 0, 0}; // velocity in ECI coordinates (km/s)
        double jd;               // julian date for time
        double el;               // elevation (degrees)
        double az;               // azimuth   (degrees)
        double range;            // range (km)

        Observer(Location location, Satellite satellite);
    };
}
#endif