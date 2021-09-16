#include <math.h>
#include <vector>
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
    CHECK( recef[0] == Approx(satellite.recef_[0]).margin(0.01).epsilon(1e-12) );
    CHECK( recef[1] == Approx(satellite.recef_[1]).margin(0.01).epsilon(1e-12) );
    CHECK( recef[2] == Approx(satellite.recef_[2]).margin(0.01).epsilon(1e-12) );
}

TEST_CASE("ComputeTeme2ecef(), skyfield test case", "[satellite][rotations]")
{
    int err;
    double jd = 2453101.827406783;
    double rteme[] = {5094.18016210, 6127.64465950, 6380.34453270};
    double recef[3];
    double recef_expected[] = {-1033.47503136, 7901.30558557, 6380.3445327};
    err = passpredict::ComputeTeme2Ecef(jd, rteme, recef);
    if ( err ) throw 1;
    CHECK( recef[0] == Approx(recef_expected[0]).margin(0.01).epsilon(1e-12) );
    CHECK( recef[1] == Approx(recef_expected[1]).margin(0.01).epsilon(1e-12) );
    CHECK( recef[2] == Approx(recef_expected[2]).margin(0.01).epsilon(1e-12) );
}

TEST_CASE("ComputeEcef2Sez(), Curtis Eg 5.9, p.271", "[satellite][rotations]")
{
    // Note that in Curtis book, the rho vector is oriented {i,j,k} => {east, north, zenith}
    // While algorithm follows Vallado orientation {i,j,k} => {south, east, zenith}
    std::vector<double> rho = {-359.0, -6.342, -466.9};
    std::vector<double> rsez;
    double lat=-40.0, lon=110.0;
    // std::vector<double> rsez_expected_curtis = {339.5, -282.6, 389.6};
    std::vector<double> rsez_expected = {282.6, 339.5, 389.6};  // for {south, east, zenith}
    rsez = passpredict::ComputeEcef2Sez(rho, lon, lat);
    CHECK( rsez[0] == Approx(rsez_expected[0]).margin(0.1).epsilon(1e-12) );
    CHECK( rsez[1] == Approx(rsez_expected[1]).margin(0.1).epsilon(1e-12) );
    CHECK( rsez[2] == Approx(rsez_expected[2]).margin(0.1).epsilon(1e-12) );
}

TEST_CASE("ComputeSez2Razel(), Curtis Eg 5.9, p.271", "[satellite][rotations]")
{
    std::vector<double> rsez = {339.5, -282.6, 389.6};
    double az, el, range;
    passpredict::ComputeSez2Razel(rsez, range, az, el);
    CHECK( range == Approx(589.0).margin(0.1).epsilon(1e-12) );
    CHECK( az == Approx(129.8).margin(0.1).epsilon(1e-12) );
    CHECK( el == Approx(41.41).margin(0.01).epsilon(1e-12) );
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

TEST_CASE("Satellite ComputeAltitude()", "[satellite]")
{
    // Vallado, Eg 3-3, p 173
    passpredict::Satellite satellite;
    satellite.recef_[0] = 6524.834;
    satellite.recef_[1] = 6862.875;
    satellite.recef_[2] = 6448.296;
    satellite.ComputeAltitude();
    REQUIRE(satellite.alt_ == Approx(5085.22));
}