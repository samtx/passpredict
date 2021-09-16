#include <array>

#include "catch.hpp"
#include "passpredict.h"


TEST_CASE( "FindAOS", "[passpredict]" ) {
    char tle1[] = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    char tle2[] = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    passpredict::Orbit orbit (tle1, tle2, wgs72);
    passpredict::Satellite satellite (orbit);
    double lat = 32.1234;
    double lon = -97.9876;
    double h = 512;
    passpredict::Location location(lat, lon, h);
    double jd0 = 2451722.5;
    double jdmax = jd0 + 14;
    double aos;
    aos = passpredict::FindAOS(jd0, jdmax, location, satellite);
    REQUIRE( aos > 0 );
    CHECK( aos > jd0 );
    CHECK( aos < jdmax );
}

