#include "catch.hpp"
#include "passpredict.h"

TEST_CASE( "Location constructor with coordinates", "[location]" ) {
    double lat = 32.1234;
    double lon = -97.9876;
    double h = 512;
    passpredict::Location location(lat, lon, h);
    REQUIRE( location.lat_ == lat );
    REQUIRE( location.lon_ == lon );
    REQUIRE( location.h_ == h );
}

TEST_CASE( "Location.site_ECEF() 39.007,-104.883,2187.0", "[location]" ) {
    // Vallado, Eg 7-1, p.431
    double lat = 39.007;
    double lon = -104.883;
    double h = 2187.0;
    passpredict::Location location(lat, lon, h);
    REQUIRE( location.recef_[0] == Approx(-1275.1219).margin(1e-4) );
    REQUIRE( location.recef_[1] == Approx(-4797.9890).margin(1e-4) );
    REQUIRE( location.recef_[2] == Approx(3994.2975).margin(1e-4) );
}

TEST_CASE( "Location.site_ECEF() 42.38,-71.13,24.0", "[location]" ) {
    double lat = 42.38;
    double lon = -71.13;
    double h = 24.0;
    passpredict::Location location(lat, lon, h);
    REQUIRE( location.recef_[0] == Approx(1526.122).margin(1e-2) );
    REQUIRE( location.recef_[1] == Approx(-4465.064).margin(1e-2) );
    REQUIRE( location.recef_[2] == Approx(4276.894).margin(1e-2) );
}