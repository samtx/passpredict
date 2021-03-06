# test rotations.py
import datetime

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

from .. import transform
from .. import sidereal
from .. import nutation
from passpredict import constants
from passpredict import timefn
from passpredict import topocentric
from passpredict.utils import epoch_from_tle_datetime
from passpredict.constants import DEG2RAD, RAD2DEG


def test_ecef2sez():
    """
    Vallado, Eg. 11-6, p.912
    """
    phi = 42.38  # latitude, deg
    lmda = -71.13  # longitude, deg
    # lmda = 136.2944
    h = 24  # height, m
    rsat = np.array([885.7296, -4389.3856, 5070.1765])
    rsite = topocentric.site_ECEF2(phi, lmda, h)
    rhoECEF = rsat - rsite
    print(rhoECEF)
    rSEZ = transform.ecef2sez(rhoECEF, phi, lmda)
    rSEZ_true = np.array([-773.8654, -581.4980, 328.8145])
    np.set_printoptions(precision=8)
    for i in [0, 1, 2]:
        assert_almost_equal(rSEZ[i], rSEZ_true[i], decimal=0, verbose=True)


@pytest.mark.parametrize(
    'jd, theta_expected, decimal',
    [
        (2448855.009722222, 152.578787810, 9),  # Vallado, Eg. 3-5
        (2453101.827406783, 312.8098943, 7)     # Vallado, Eg. 3-15, p.230

    ],
    ids=[
        'Vallado Eg. 3-5',
        'Vallado Eg. 3-15',
    ]
)
def test_theta_GMST1982(jd, theta_expected, decimal):
    """Compute the Greenwich Mean Sidereal Time

    References:
        Vallado, Eg 3-5, p.188
    """
    # dt = datetime.datetime(1992, 8, 20, 12, 14, 0)  # Aug 20, 1992, 12:14 PM UT1
    # jd = timefn.julian_date(dt)  # jd = 2448855.009722222
    theta, thetadt = sidereal.theta_GMST1982(jd)
    theta *= constants.RAD2DEG  # convert from radians to degrees
    assert_almost_equal(theta, theta_expected, decimal=decimal)  # degrees


def test_IJK2SEZ():
    """
    Curtis, Eg. 5.9, p.270
    """
    from math import sin, cos, radians, asin, acos, degrees

    # Satellite position
    r = np.array([-2032.4, 4591.2, -4544.8])  # km, geocentric equatorial pos.
    # Location position
    H = 0  # elevation, sea level
    phi = -40  # deg latitude
    theta = 110  # deg, local sidereal time
    R_obs = np.array([-1673.0, 4598.0, -4078.0])  # km, from Eq. 5.56

    rho = r - R_obs  # relative position vector
    phi_rad = radians(phi)
    theta_rad = radians(theta)

    # rotation matrix, from Eq. 5.62a
    Q = np.array(
        [
            [-sin(theta_rad), cos(theta_rad), 0.0],
            [
                -sin(phi_rad) * cos(theta_rad),
                -sin(phi_rad) * sin(theta_rad),
                cos(phi_rad),
            ],
            [
                cos(phi_rad) * cos(theta_rad),
                cos(phi_rad) * sin(theta_rad),
                sin(phi_rad),
            ],
        ]
    )

    rho2 = np.dot(Q, rho)

    rho2_hat = (1 / np.linalg.norm(rho)) * rho2

    a = degrees(asin(rho2_hat[2]))
    A = degrees(acos(rho2_hat[1] / cos(radians(a))))

    A_true = 129.8  # deg, azimuth
    a_true = 41.41  # deg, angular elevation
    assert_almost_equal(A, A_true, decimal=1)
    assert_almost_equal(a, a_true, decimal=1)




@pytest.mark.skip()
def test_appendix_c_conversion_from_TEME_to_ITRF_UTC1():
    """Test TEME to ITRF conversion

    References:
        Vallado et al., Revision 2
        Rhodes, Skyfield library, test_earth_satellites.py
    """
    seconds_per_day = 24.0 * 60.0 * 60.0
    rTEME = np.array([5094.18016210, 6127.64465950, 6380.34453270])
    vTEME = np.array([-4.746131487, 0.785818041, 5.531931288])
    vTEME = vTEME * seconds_per_day  # km/s to km/day

    # Apr 6, 2004,  07:51:28.386 UTC
    deltaUTC1 = -439961 # microseconds
    ms = int(1e6) + 386000# + deltaUTC1
    ms = 386000 - 439961 + int(1e6)
    dt = datetime.datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)

    # Polar motion
    xp = -0.140682  # arcseconds
    yp = 0.333309  # arcseconds
    xp *= constants.ASEC2RAD
    yp *= constants.ASEC2RAD
    rITRF, vITRF = rotations.TEME_to_ITRF(jd, rTEME, vTEME, xp, yp)

    print(rITRF)
    assert_almost_equal(rITRF[0], -1033.47938300, decimal=4)
    assert_almost_equal(rITRF[1], 7901.29527540, decimal=4)
    assert_almost_equal(rITRF[2], 6380.35659580, decimal=4)

    vITRF /= seconds_per_day  # km/day to km/s
    print(vITRF)
    assert_almost_equal(vITRF[0], -3.225636520, decimal=6)
    assert_almost_equal(vITRF[1], -2.872451450, decimal=6)
    assert_almost_equal(vITRF[2], 5.531924446, decimal=6)


# def test_appendix_c_conversion_from_TEME_to_ITRF_with_teme2ecef():
#     """Test TEME to ITRF conversion

#     References:
#         Vallado et al., Revision 2
#         Rhodes, Skyfield library, test_earth_satellites.py
#     """
#     rTEME = np.array([[5094.18016210], [6127.64465950], [6380.34453270]])
    
#     # Apr 6, 2004,  07:51:28.386 UTC
#     seconds = 28.386
#     deltaUTC1 = -0.439961 # seconds
#     s, us = np.divmod(seconds + deltaUTC1, 1)
#     us *= 1e6  # microseconds
#     dt = datetime.datetime(2004, 4, 6, 7, 51, int(s), int(us))
#     jd = timefn.julian_date(dt)

#     # Polar motion
#     xp = -0.140682  # arcseconds
#     yp = 0.333309  # arcseconds
#     xp *= constants.ASEC2RAD
#     yp *= constants.ASEC2RAD
#     rITRF = transform.teme2ecef(rTEME, jd, xp, yp)

#     print(rITRF)
#     assert_almost_equal(rITRF[0], -1033.47938300, decimal=4)
#     assert_almost_equal(rITRF[1], 7901.29527540, decimal=4)
#     assert_almost_equal(rITRF[2], 6380.35659580, decimal=4)


if __name__ == "__main__":
    import pytest
    pytest.main(['-v', __file__])