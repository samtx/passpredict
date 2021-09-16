#include <math.h>

#include "location.h"
#include "passpredict.h"

namespace passpredict {

Location::Location(double lat, double lon, double h) : lat_(lat), lon_(lon), h_(h)
{
    Location::SiteEcef();
};

// get ECEF position vector
void Location::SiteEcef()
{
    //  Vallado, Eq. 3-7
    double phi_gd_rad, h_ellp_km, C, S, r_delta, r_K, sin_phi_gd, lmda_rad;
    phi_gd_rad = lat_ * PASSPREDICT_DEG2RAD;
    sin_phi_gd = std::sin(phi_gd_rad);
    h_ellp_km = h_ * 0.001; // convert height to km
    C = PASSPREDICT_R_EARTH / std::sqrt(1 - PASSPREDICT_e2_EARTH * std::pow(sin_phi_gd, 2));
    S = C * (1 - PASSPREDICT_e2_EARTH);
    r_delta = (C + h_ellp_km) * std::cos(phi_gd_rad);
    r_K = (S + h_ellp_km) * sin_phi_gd;

    // Vallado, Alg 51, p.430
    lmda_rad = lon_ * PASSPREDICT_DEG2RAD;
    recef_[0] = r_delta * std::cos(lmda_rad);
    recef_[1] = r_delta * std::sin(lmda_rad);
    recef_[2] = r_K;
};

} // namespace passpredict
