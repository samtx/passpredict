#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include "catch.hpp"
#include "SGP4.h"
#include "passpredict.h"

// #include <catch2/catch_test_macros.hpp>
// #include <catch2/generators/catch_generators.hpp>
// #include <catch2/generators/catch_generators_adapters.hpp>

// Custom Catch2 generator for comparing TLE and propagation output

// Ref: https://github.com/catchorg/Catch2/blob/devel/examples/302-Gen-Table.cpp
// Ref: https://github.com/catchorg/Catch2/blob/devel/examples/300-Gen-OwnGenerator.cpp

// class Sgp4SatGenerator : public Catch::Generators::IGenerator<int> {
//     int m_satid;

// public:
//     Sgp4SatGenerator(int satid): m_satid(satid) {

//     }
// }

// Loop over TLEs in SGP4-VER.TLE file and compare with output in tcppverDec2015.out
TEST_CASE( "Satellite SGP4 Propagation, satid 5", "[satellite][sgp4]" ) {
    // From SGP4-VER.TLE and tcppverDec2015.out
    char tle1[] = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    char tle2[] = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    double tsince[] = {0, 360, 720, 1080, 1440, 1800, 2160, 2520, 2880, 3240, 3600, 3960, 4320};
    double r[][3] = {
        {  7022.46529266, -1400.08296755,     0.03995155},
        { -7154.03120202, -3783.17682504, -3536.19412294},
        { -7134.59340119,  6531.68641334,  3260.27186483},
        {  5568.53901181,  4492.06992591,  3863.87641983},
        {  -938.55923943, -6268.18748831, -4294.02924751},
        { -9680.56121728,  2802.47771354,   124.10688038},
        {   190.19796988,  7746.96653614,  5110.00675412},
        {  5579.55640116, -3995.61396789, -1518.82108966},
        { -8650.73082219, -1914.93811525, -3007.03603443},
        { -5429.79204164,  7574.36493792,  3747.39305236},
        {  6759.04583722,  2001.58198220,  2783.55192533},
        { -3791.44531559, -5712.95617894, -4533.48630714},
        { -9060.47373569,  4658.70952502,   813.68673153}
    };
    double v[][3] = {
        {  1.893841015,  6.405893759,  4.534807250 },
        {  4.741887409, -4.151817765, -2.093935425 },
        { -4.113793027, -2.911922039, -2.557327851 },
        { -4.209106476,  5.159719888,  2.744852980 },
        {  7.536105209, -0.427127707,  0.989878080 },
        { -0.905874102, -4.659467970, -3.227347517 },
        { -6.112325142,  1.527008184, -0.139152358 },
        {  4.767927483,  5.123185301,  4.276837355 },
        {  3.067165127, -4.828384068, -2.515322836 },
        { -4.999442110, -1.800561422, -2.229392830 },
        { -2.180993947,  6.402085603,  3.644723952 },
        {  6.668817493, -2.516382327, -0.082384354 },
        { -2.232832783, -4.110453490, -3.157345433 }
    };

    passpredict::Orbit orbit (tle1, tle2, wgs72);
    passpredict::Satellite satellite (orbit);
    int i = GENERATE(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
    satellite.tsince = tsince[i];
    satellite.sgp4();
    REQUIRE( satellite.rteme[0] == Approx(r[i][0]) );
    REQUIRE( satellite.rteme[1] == Approx(r[i][1]) );
    REQUIRE( satellite.rteme[2] == Approx(r[i][2]) );
    REQUIRE( satellite.vteme[0] == Approx(v[i][0]) );
    REQUIRE( satellite.vteme[1] == Approx(v[i][1]) );
    REQUIRE( satellite.vteme[2] == Approx(v[i][2]) );
}
