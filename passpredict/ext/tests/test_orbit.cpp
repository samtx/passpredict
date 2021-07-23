#include "catch.hpp"
#include "passpredict.h"

TEST_CASE( "Orbit constructor with TLE strings", "[orbit]" ) {
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    passpredict::Orbit orbit(tle1, tle2);
    REQUIRE( orbit.satrec.jdsatepoch == 2459415.5 );
    REQUIRE( orbit.satrec.jdsatepochF == Approx(0.469801) );
    REQUIRE( orbit.satrec.bstar == 4.2487e-5 );
    REQUIRE( orbit.satrec.inclo / PASSPREDICT_DEG2RAD == 51.6426 );
    REQUIRE( orbit.satrec.nodeo / PASSPREDICT_DEG2RAD == 178.1369 );
    REQUIRE( orbit.satrec.ecco == 0.0001717 );
    REQUIRE( orbit.satrec.argpo / PASSPREDICT_DEG2RAD == 174.7410 );
    REQUIRE( orbit.satrec.mo / PASSPREDICT_DEG2RAD == 330.7918 );
    REQUIRE( orbit.satrec.no_kozai * (1440.0 / PASSPREDICT_2PI) == Approx(15.4883) );
}
