#include <math.h>
#include <vector>
#include "catch.hpp"
extern "C"
{
#include "sofa.h"
}
#include "satellite.h"
#include "passpredict.h"

TEST_CASE("Satellite Teme2Ecef(), skyfield test case", "[satellite][rotations]")
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
    std::vector<double> rsez = {282.6, 339.5, 389.6};
    double az, el, range;
    passpredict::ComputeSez2Razel(rsez, range, az, el);
    CHECK( range == Approx(589.0).margin(0.1).epsilon(1e-12) );
    CHECK( az == Approx(129.8).margin(0.1).epsilon(1e-12) );
    CHECK( el == Approx(41.41).margin(0.01).epsilon(1e-12) );
}

TEST_CASE("Rho ECEF vector to Razel coords, Curtis Eg 5.9, p.271", "[satellite][rotations]")
{
    // Note that in Curtis book, the horizon vector is oriented {i,j,k} => {east, north, zenith}
    // While algorithm follows Vallado orientation {i,j,k} => {south, east, zenith}
    std::vector<double> rho = {-359.0, -6.342, -466.9};
    std::vector<double> rsez;
    double az, el, range;
    double lat=-40.0, lon=110.0;
    // std::vector<double> rsez_expected_curtis = {339.5, -282.6, 389.6};
    double range_expected = 589.0;
    double az_expected = 129.8;
    double el_expected = 41.41;
    rsez = passpredict::ComputeEcef2Sez(rho, lon, lat);
    std::vector<double> rsez_expected = {282.6, 339.5, 389.6};  // for {south, east, zenith}
    CHECK( rsez[0] == Approx(rsez_expected[0]).margin(0.1).epsilon(1e-12) );
    CHECK( rsez[1] == Approx(rsez_expected[1]).margin(0.1).epsilon(1e-12) );
    CHECK( rsez[2] == Approx(rsez_expected[2]).margin(0.1).epsilon(1e-12) );
    passpredict::ComputeSez2Razel(rsez, range, az, el);
    CHECK( range == Approx(range_expected).margin(0.1).epsilon(1e-12) );
    CHECK( az == Approx(az_expected).margin(0.1).epsilon(1e-12) );
    CHECK( el == Approx(el_expected).margin(0.01).epsilon(1e-12) );
}

TEST_CASE("Rho ECEF vector to Razel coords, Vallado Eg 11-6, p.913", "[satellite][rotations]")
{
    // View data table 11-4 for results for Apr 2, 1997 1:08:00 UTC
    passpredict::Location location;
    std::vector<double> satellite_ecef = {885.7296, -4389.3856, 5070.1765};
    std::vector<double> rho_sez(3, 0.0);
    std::vector<double> rho_ecef(3, 0.0);
    int i;
    double az, el, range;
    double lat=42.38, lon=-71.13, h=24;
    location = passpredict::Location(lat, lon, h);
    for (i=0; i<3; i++)
        rho_ecef[i] = satellite_ecef[i] - location.recef_[i];
    rho_sez = passpredict::ComputeEcef2Sez(rho_ecef, lon, lat);
    CHECK( rho_sez[0] == Approx(-773.8654).margin(0.0001).epsilon(1e-12) );
    CHECK( rho_sez[1] == Approx(-581.4980).margin(0.0001).epsilon(1e-12) );
    CHECK( rho_sez[2] == Approx(328.8145).margin(0.0001).epsilon(1e-12) );
    passpredict::ComputeSez2Razel(rho_sez, range, az, el);
    CHECK( range == Approx(1022.3143).margin(0.0001).epsilon(1e-12) );
    CHECK( az == Approx(323.0780).margin(0.0001).epsilon(1e-12) );
    CHECK( el == Approx(18.7619).margin(0.0001).epsilon(1e-12) );
}

TEST_CASE("ComputeAzimuth() {282.6, 339.5, 389.6}", "[satellite][rotations]"){
    // Curtis example
    double rsez[3] = {282.6, 339.5, 389.6};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(129.8).margin(0.1).epsilon(1e-12) );
}

TEST_CASE("ComputeAzimuth() {-773.8654, -581.4980, 328.8145}", "[satellite][rotations]"){
    // Vallado example 11-6
    double rsez[3] = {-773.8654, -581.4980, 328.8145};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(323.0780).margin(0.0001).epsilon(1e-12) );
}

TEST_CASE("ComputeAzimuth() 45 deg NE", "[satellite][rotations]"){
    //
    double rsez[3] = {-1000, 1000, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(45.00) );
}

TEST_CASE("ComputeAzimuth() 135 deg SE", "[satellite][rotations]"){
    //
    double rsez[3] = {1000, 1000, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(135.00) );
}

TEST_CASE("ComputeAzimuth() 225 deg SW", "[satellite][rotations]"){
    //
    double rsez[3] = {1000, -1000, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(225.00) );
}

TEST_CASE("ComputeAzimuth() 315 deg SW", "[satellite][rotations]"){
    //
    double rsez[3] = {-1000, -1000, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(315.00) );
}

TEST_CASE("ComputeAzimuth() 90 deg east", "[satellite][rotations]"){
    //
    double rsez[3] = {0, 1000, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(90.00) );
}

TEST_CASE("ComputeAzimuth() 180 deg south", "[satellite][rotations]"){
    //
    double rsez[3] = {1000, 0, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(180.00) );
}


TEST_CASE("ComputeAzimuth() 270 deg west", "[satellite][rotations]"){
    //
    double rsez[3] = {0, -1000, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(270.00) );
}

TEST_CASE("ComputeAzimuth() 0 deg north", "[satellite][rotations]"){
    //
    double rsez[3] = {-1000, 0, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(0.00) );
}

TEST_CASE("ComputeAzimuth() 0.05 deg northeast", "[satellite][rotations]"){
    //
    double rsez[3] = {-1000, 1, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(0.05729).margin(0.0001).epsilon(1e-12) );
}

TEST_CASE("ComputeAzimuth() 359.94 deg northwest", "[satellite][rotations]"){
    //
    double rsez[3] = {-1000, -1, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(359.942704).margin(0.0001).epsilon(1e-12) );
}

TEST_CASE("ComputeAzimuth() 179.9427 deg southeast", "[satellite][rotations]"){
    //
    double rsez[3] = {1000, 1, 1000};
    double az;
    az = passpredict::ComputeAzimuth(rsez[0], rsez[1]);
    CHECK( az == Approx(179.9427).margin(0.0001).epsilon(1e-12) );
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