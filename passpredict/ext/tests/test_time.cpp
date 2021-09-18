#include "catch.hpp"
extern "C"
{
#include "sofa.h"
}
#include "satellite.h"
#include "passpredict.h"

TEST_CASE("Time Utc2tt(), SOFA example, p.3", "[satellite][time]")
{
    int year = 2010, month = 7, day = 24, hour = 11, minute = 18;
    int iy, im, id, ihmsf[4];
    double second = 7.318;
    double jd1, jd2, tt1, tt2;
    int j;
    j = iauDtf2d("UTC", year, month, day, hour, minute, second, &jd1, &jd2);
    REQUIRE( !j );
    REQUIRE( jd1 == Approx(2455401.5000) );
    REQUIRE( jd2 == Approx(0.470918032407) );
    // double jd1 = 2455401.5000, jd2 = 0.470918032407;
    j = passpredict::Utc2tt(jd1, jd2, tt1, tt2);
    REQUIRE( !j );
    CHECK( tt1 == Approx(2455401.5) );
    CHECK( tt2 == Approx(0.471684050926) );
    j = iauD2dtf("tt", 3, tt1, tt2, &iy, &im, &id, ihmsf);
    REQUIRE( !j );
    CHECK( iy == 2010 );
    CHECK( im == 7 );
    CHECK( id == 24 );
    CHECK( ihmsf[0] == 11 );
    CHECK( ihmsf[1] == 19 );
    CHECK( ihmsf[2] == 13 );
    CHECK( ihmsf[3] == 502 );
}