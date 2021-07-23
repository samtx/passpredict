#include "catch.hpp"
#include "passpredict.h"

TEST_CASE( "Location constructor with coordinates", "[location]" ) {
    double lat = 32.1234;
    double lon = -97.9876;
    double h = 512;
    passpredict::Location location(lat, lon, h);
    REQUIRE( location.lat == lat );
    REQUIRE( location.lon == lon );
    REQUIRE( location.h == h );
}