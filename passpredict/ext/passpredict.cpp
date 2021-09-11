#include <array>

#include "satellite.h"
#include "location.h"

#include "passpredict.h"


namespace passpredict {

// double ComputeElevationAngle(double jd, Location location, Satellite satellite){
//     // Compute elevation angle from location to satellite at time jd
//     double el;
//     int i;

//     std::array<double, 3> rho;
//     std::array<double, 3> rsez;

//     // Propagate satellite to time jd_
//     sat_ptr_->PropagateJd(jd);

//     // get vector rho from location to satellite
//     for (i=0; i<3; i++)
//         rho[i] = sat_ptr_->recef_[i] - loc_ptr_->recef_[i];

//     // Rotate rho ECEF vector to topographic SEZ vector
//     Observer::Ecef2Sez(rho, rsez);

//     // Compute elevation angle, range, and azimuth
//     Observer::Sez2Razel(rsez);

//     return el;
// }



} // passpredict