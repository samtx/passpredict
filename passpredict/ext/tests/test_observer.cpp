#include "catch.hpp"
#include "passpredict.h"

//     phi = 42.38  # latitude, deg
//     lmda = -71.13  # longitude, deg
//     # lmda = 136.2944
//     h = 24  # height, m
//     rsat = np.array([885.7296, -4389.3856, 5070.1765])
//     rsite = topocentric.site_ECEF2(phi, lmda, h)
//     rhoECEF = rsat - rsite
//     print(rhoECEF)
//     rSEZ = rotations.ecef2sez(rhoECEF, phi, lmda)
//     rSEZ_true = np.array([-773.8654, -581.4980, 328.8145])
//     np.set_printoptions(precision=8)
//     assert_allclose(rSEZ, rSEZ_true)


TEST_CASE( "Observer Ecef2Sez", "[observer]" ) {
    double rsat[3] = {885.7296, -4389.3856, 5070.1765}; // ecef
    double rsez[3], rho[3];
    int i;
    passpredict::Location location(42.38, -71.13, 24);
    passpredict::Satellite satellite;
    passpredict::Observer observer(location, satellite);
    for (i=0; i<3; i++)
        rho[i] = rsat[i] - location.recef_[i];
    observer.Ecef2Sez(rho, rsez);
    CHECK( rsez[0] == Approx(-773.8654) );
    CHECK( rsez[1] == Approx(-581.4980) );
    CHECK( rsez[2] == Approx(328.8145) );
}

