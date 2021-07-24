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
        satellite.rteme[i] = rteme[i];
    satellite.jd = 2453101.827406783;
    double recef[] = {-1033.47503136, 7901.30558557, 6380.3445327};
    satellite.teme2ecef();
    CHECK(recef[0] == Approx(satellite.recef[0]));
    CHECK(recef[1] == Approx(satellite.recef[1]));
    CHECK(recef[2] == Approx(satellite.recef[2]));

}