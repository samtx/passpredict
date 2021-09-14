#include <array>
#include <vector>
#include <cmath>

#include "satellite.h"
#include "location.h"
#include "passpredict.h"


namespace passpredict {

double ComputeElevationAngle(double jd, Location location, Satellite satellite){
    // Compute elevation angle from location to satellite at time jd
    double el, range;
    int i;

    std::vector<double> rho(3, 0.0), rsez(3, 0.0), recef(3, 0.0);

    // Propagate satellite to time jd_
    recef = PropagateSatelliteJd(jd, satellite);

    // get vector rho from location to satellite
    for (i=0; i<3; i++)
        rho[i] = recef[i] - location.recef_[i];

    // Rotate rho ECEF vector to topographic SEZ vector
    rsez = ComputeEcef2Sez(rho, location.lon_, location.lat_);

    // Compute elevation angle from SEZ vector
    range = std::sqrt(rsez[0]*rsez[0] + rsez[1]*rsez[1] + rsez[2]*rsez[2]);
    el = std::asin(rsez[2] / range) * PASSPREDICT_RAD2DEG;

    return el;
}


double FindAOS(double jd0, double jdmax, Location location, Satellite satellite){
    // find next AOS time
    double jd, jdstep, el, elold, jdold, el1, el2, jd1, jd2, tmp;
    double dt_second = 1 / 86400.0;
    int i;

    jd = jd0;
    el = -90.0;

    if (jdmax <= 0.0)
        // error
        return -1;
    // coarse time steps, 20 seconds
    while ((jd < jdmax) && (el < 0)) {
        elold = el;
        jdold = jd;
        jd += 20.0/86400.0;  // add 20 seconds
        el = ComputeElevationAngle(jd, location, satellite);
    }
    if ((el < 0) || (elold > 0)) {
        // error
        return -1;
    }
    // fine time steps, 1 second
    tmp = jdold;
    jd2 = jd;
    jd = tmp;
    while ((jd <= jd2) && (el < 0)) {
        elold = el;
        jdold = jd;
        jd += 1.0/86400.0;  // add 1 second
        el = ComputeElevationAngle(jd, location, satellite);
    }
    if ((elold < 0) && (el >= 0)) {
        // time of AOS found
        return jd;
    }
    // error
    return -1;
}


} // passpredict