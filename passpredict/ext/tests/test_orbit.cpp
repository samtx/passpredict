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

TEST_CASE( "Orbit constructor with OMM struct", "[orbit]" ) {
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    passpredict::Orbit orbit_tle(tle1, tle2, wgs72);
    passpredict::Omm omm;
    omm.jdsatepoch = 2459415.5;
    omm.jdsatepochF = 0.469801;
    omm.bstar = 4.2487e-5;
    omm.inclo = 51.6426;
    omm.ecco = 0.0001717;
    omm.argpo = 174.7410;
    omm.mo = 330.7918;
    omm.nodeo = 178.1369;
    omm.no_kozai = 15.4883;
    omm.elnum = 9993;
    omm.revnum = 293750;
    strcpy(omm.satnum, "25544");
    omm.classification = 'U';
    omm.ephtype = 0;
    passpredict::Orbit orbit_omm(omm);

    CHECK( orbit_tle.satrec_.jdsatepoch == orbit_omm.satrec_.jdsatepoch );
    CHECK( orbit_tle.satrec_.jdsatepochF == Approx(orbit_omm.satrec_.jdsatepochF).margin(1e-6).epsilon(1e-12) );
    CHECK( orbit_tle.satrec_.bstar == orbit_omm.satrec_.bstar );
    CHECK( orbit_tle.satrec_.inclo == orbit_omm.satrec_.inclo );
    CHECK( orbit_tle.satrec_.nodeo == orbit_omm.satrec_.nodeo );
    CHECK( orbit_tle.satrec_.ecco == orbit_omm.satrec_.ecco );
    CHECK( orbit_tle.satrec_.argpo == orbit_omm.satrec_.argpo );
    CHECK( orbit_tle.satrec_.mo == orbit_omm.satrec_.mo );
    CHECK( orbit_tle.satrec_.no_kozai == Approx(orbit_omm.satrec_.no_kozai).margin(1e-4).epsilon(1e-12) );
}


TEST_CASE("Orbit constructor with twoline2rv from sgp4", "[sgp4][satellite][orbit]"){
    char tle1[] = "1 25544U 98067A   21201.46980141  .00001879  00000-0  42487-4 0  9993";
    char tle2[] = "2 25544  51.6426 178.1369 0001717 174.7410 330.7918 15.48826828293750";
    char tle1_copy[strlen(tle1) + 1];
    char tle2_copy[strlen(tle2) + 1];
    strcpy(tle1_copy, tle1);
    strcpy(tle2_copy, tle2);
    REQUIRE( !strcmp(tle1, tle1_copy) );
    REQUIRE( !strcmp(tle2, tle2_copy) );
    passpredict::Orbit orbit (tle1, tle2);
    // char typerun = ' ';
    // char typeinput = ' ';
    // char opsmode = 'i';
    gravconsttype whichconst = wgs84;
    elsetrec satrec;
    double dummy;
    SGP4Funcs::twoline2rv(tle1_copy, tle2_copy, ' ', ' ', 'i',
        whichconst, dummy, dummy, dummy, satrec
    );
    CHECK( orbit.satrec_.nodeo == satrec.nodeo );
    CHECK( orbit.satrec_.bstar == satrec.bstar );
    CHECK( orbit.satrec_.argpo == satrec.argpo );
    CHECK( orbit.satrec_.ecco == satrec.ecco );
    CHECK( orbit.satrec_.inclo == satrec.inclo );
    CHECK( orbit.satrec_.mo == satrec.mo );
    CHECK( orbit.satrec_.no_kozai == satrec.no_kozai );
    CHECK( orbit.satrec_.jdsatepoch == satrec.jdsatepoch );
    CHECK( orbit.satrec_.jdsatepochF == satrec.jdsatepochF );
}