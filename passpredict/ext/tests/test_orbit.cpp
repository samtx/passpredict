#include "catch.hpp"
#include "SGP4.h"
#include "passpredict.h"

TEST_CASE( "Orbit constructor with TLE strings", "[orbit]" ) {
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    passpredict::Orbit orbit(tle1, tle2);
    CHECK( orbit.satrec_.jdsatepoch == 2459415.5 );
    CHECK( orbit.satrec_.jdsatepochF == Approx(0.469801) );
    CHECK( orbit.satrec_.bstar == 4.2487e-5 );
    CHECK( orbit.satrec_.inclo / PASSPREDICT_DEG2RAD == 51.6426 );
    CHECK( orbit.satrec_.nodeo / PASSPREDICT_DEG2RAD == 178.1369 );
    CHECK( orbit.satrec_.ecco == 0.0001717 );
    CHECK( orbit.satrec_.argpo / PASSPREDICT_DEG2RAD == 174.7410 );
    CHECK( orbit.satrec_.mo / PASSPREDICT_DEG2RAD == 330.7918 );
    CHECK( orbit.satrec_.no_kozai * (1440.0 / PASSPREDICT_2PI) == Approx(15.4883) );
}


TEST_CASE( "Orbit constructor with TLE strings and wgs72", "[orbit]" ) {
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    passpredict::Orbit orbit(tle1, tle2, wgs72);
    CHECK( orbit.satrec_.jdsatepoch == 2459415.5 );
    CHECK( orbit.satrec_.jdsatepochF == Approx(0.469801) );
    CHECK( orbit.satrec_.bstar == 4.2487e-5 );
    CHECK( orbit.satrec_.inclo / PASSPREDICT_DEG2RAD == 51.6426 );
    CHECK( orbit.satrec_.nodeo / PASSPREDICT_DEG2RAD == 178.1369 );
    CHECK( orbit.satrec_.ecco == 0.0001717 );
    CHECK( orbit.satrec_.argpo / PASSPREDICT_DEG2RAD == 174.7410 );
    CHECK( orbit.satrec_.mo / PASSPREDICT_DEG2RAD == 330.7918 );
    CHECK( orbit.satrec_.no_kozai * (1440.0 / PASSPREDICT_2PI) == Approx(15.4883) );
}