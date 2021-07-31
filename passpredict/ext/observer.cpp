#include <math.h>
#include "observer.h"
#include "passpredict.h"

namespace passpredict {

Observer::Observer(Location location, Satellite satellite) 
    : location_(location), satellite_(satellite) {};

double Observer::ComputeElevationAngle(){
    // Compute observation elevation angle at time jd_
    int i;
    double rho[3], rsez[3];
    
    // Propagate satellite to time jd_
    satellite_.PropagateJd(jd_);

    // get vector rho from location to satellite
    for (i=0; i<3; i++)
        rho[i] = satellite_.recef_[i] - location_.recef_[i];

    // Rotate rho ECEF vector to topographic SEZ vector 
    Observer::Ecef2Sez(rho, rsez);

    // Compute elevation angle, range, and azimuth
    Observer::Sez2Razel(rsez);

    return el_;

};

// void Observer::FindAos(double t0, double tmax) {
//     // find r,az,el at t0
// };

// void Observer::FindLos() {};

void Observer::Ecef2Sez(double recef[3], double rsez[3]) {
    // Rotate vector from ECEF to SEZ frame based on 
    // location geodetic coordinates
    double phi_rad, lmda_rad, ang1, cosang1, sinang1;
    double cosang2, sinang2;
    phi_rad = location_.lat_ * PASSPREDICT_DEG2RAD;
    lmda_rad = location_.lon_ * PASSPREDICT_DEG2RAD;
    ang1 = (90 - location_.lat_) * PASSPREDICT_DEG2RAD;
    cosang1 = std::cos(ang1);
    sinang1 = std::sin(ang1);
    cosang2 = std::cos(lmda_rad);
    sinang2 = std::sin(lmda_rad);
    rsez[0] = cosang1*cosang2*recef[0] + cosang1*sinang2*recef[1] - sinang2*recef[2];
    rsez[1] = -sinang2*recef[0] + cosang2*recef[1] + 0;
    rsez[2] = sinang1*cosang2*recef[0] + sinang1*sinang2*recef[1] + cosang1*recef[2];
};

void Observer::Sez2Razel(double rsez[3]) {
    // Get Range, Elevation, and Azimuth from SEZ vector
    range_ = std::sqrt(
        rsez[0]*rsez[0] + rsez[1]*rsez[1] + rsez[2]*rsez[2]
    );
    el_ = std::asin(rsez[2] / range_) * PASSPREDICT_RAD2DEG;
    az_ = std::atan2(rsez[0], rsez[1]) + PASSPREDICT_RAD2DEG;
    if (rsez[0] < 0 && rsez[1] < 0)
        az_ = std::fmod(az_, 360.0);
};


}; // namespace passpredict