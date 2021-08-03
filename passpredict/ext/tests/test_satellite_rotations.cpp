#include "catch.hpp"
extern "C"
{
#include "sofa.h"
}
#include "satellite.h"
#include "passpredict.h"

TEST_CASE("Satellite teme2ecef(), skyfield test case", "[satellite][rotations]")
{
    passpredict::Satellite satellite;
    double rteme[] = {5094.18016210, 6127.64465950, 6380.34453270};
    for (int i=0; i<3; i++)
        satellite.rteme_[i] = rteme[i];
    satellite.jd_ = 2453101.827406783;
    double recef[] = {-1033.47503136, 7901.30558557, 6380.3445327};
    satellite.Teme2Ecef();
    CHECK(recef[0] == Approx(satellite.recef_[0]));
    CHECK(recef[1] == Approx(satellite.recef_[1]));
    CHECK(recef[2] == Approx(satellite.recef_[2]));
}

TEST_CASE("Satellite ComputeSubpoint()", "[satellite]")
{
    // Vallado, Eg 3-3, p 173
    passpredict::Satellite satellite;
    passpredict::Subpoint subpoint;
    satellite.recef_[0] = 6524.834;
    satellite.recef_[1] = 6862.875;
    satellite.recef_[2] = 6448.296;
    subpoint = satellite.ComputeSubpoint();
    CHECK(subpoint.lat == Approx(34.352496));
    CHECK(subpoint.lon == Approx(46.4464));
    CHECK(subpoint.alt == Approx(5085.22));
}