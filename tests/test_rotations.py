# test rotations.py
import datetime

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

from passpredict import rotations
from passpredict import constants
from passpredict import timefn
from passpredict import topocentric
from passpredict.utils import epoch_from_tle_datetime
from passpredict import nutation
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
    rSEZ = rotations.ecef2sez(rhoECEF, phi, lmda)
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
    theta, thetadt = rotations.theta_GMST1982(jd)
    theta *= constants.RAD2DEG  # convert from radians to degrees
    assert_almost_equal(theta, theta_expected, decimal=decimal)  # degrees


def test_ecef2eci():
    """
    Vallado, Eg. 3-15, p. 230
    """
    rECEF = np.array([
        [-1033.4793830],
         [7901.2952754],
         [6380.3565958]
    ])
    jdt = timefn.julian_date(2004, 4, 6, 7, 51, 28.386009)
    rECI = rotations.ecef2eci(rECEF, jdt)
    assert_allclose(
        rECI,
        np.array([[5102.5096], [6123.01152], [6378.1363]]),
        atol=0.1  # need to make more accurate
    )


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


def test_thetaGMST_from_skyfield():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime.datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)
    theta, thetadt = rotations.theta_GMST1982(jd)
    assert_almost_equal(theta, 5.459562584754709, decimal=15)


@pytest.mark.parametrize(
    "dpsi, meaneps, eq_expected",
    [
        (-0.004250260*DEG2RAD, 23.43922657*DEG2RAD, -0.0038995*DEG2RAD),
        (-0.004337544*DEG2RAD, 23.43922657*DEG2RAD, -0.0039796*DEG2RAD),
    ],
    ids=[
        'epoch=00182.78495062',
        'epoch=00179.78495062',
    ]
)
def test_equinox1982_geometric_terms(dpsi, meaneps, eq_expected):
    """
    Vallado, p.234, epoch datetime = 00179.78495062
    """
    eq = rotations.equinox1982_geometric_terms(dpsi, meaneps)
    assert_almost_equal(eq, eq_expected)


@pytest.mark.xfail(reason='Need to fix this later')
def test_equinox1982_geometric_terms_3():
    """
    Vallado, Eg.3-15, p.230,
    """
    theta_GMST1982 = 312.8098943  # degrees
    eq_true =  312.8067654 - theta_GMST1982
    jd = timefn.julian_date(2006, 4, 6, 7, 51, 28.386009)
    tt = timefn.jd2jc(jd)
    nut = nutation.fk5_nutation(tt)
    eq = rotations.equinox1982(nut.dpsi, nut.eps, nut.omega_moon)
    assert_almost_equal(eq, eq_true * constants.DEG2RAD)


@pytest.mark.parametrize(
    'epoch_string, r0_expected, r1_expected, r2_expected',
    [
        ('00182.78495062', -9059.9413786, 4659.6972000, 813.9588875),  # Vallado, p. 234
        ('00179.78495062', -9059.9510799, 4659.6807556, 813.9450451),  # Vallado, p. 234
    ],
    ids=[
        'epoch=00182.78495062',
        'epoch=00179.78495062',
    ]
)
def test_teme2eci(epoch_string, r0_expected, r1_expected, r2_expected):
    """
    Test TEME -> J2000 (ECI) conversion

    Using 'of date'

    Reference:
        Vallado, p.234, Eq. 3-91
    """
    epoch = epoch_from_tle_datetime(epoch_string)
    print(f'epoch = {epoch}')
    jd = timefn.julian_date(epoch)
    print(f'jd = {jd}')
    rTEME = np.array((-9060.47373569, 4658.70952502, 813.68673153))
    rJ2000 = rotations.teme2eci(rTEME, jd)
    print(f'rJ2000 = {rJ2000}')
    assert_almost_equal(rJ2000[0], r0_expected, decimal=3)
    assert_almost_equal(rJ2000[1], r1_expected, decimal=3)
    assert_almost_equal(rJ2000[2], r2_expected, decimal=3)


@pytest.mark.skip('failing')
def test_theta_GAST1982():
    """
    Vallado, Eg. 3-15, p. 230
    """
    dut1 = -0.4399619
    jdut1 = timefn.julian_date(2004, 4, 6, 7, 51, 28.386009 + dut1)
    ttut1 = timefn.jd2jc(jdut1)
    gmst = 312.8098943 * constants.DEG2RAD
    tt = 0.0426236319
    nut = nutation.fk5_nutation.__wrapped__(tt)
    print('dpsi ',nut.dpsi*constants.RAD2DEG)
    print('eps ',nut.eps*constants.RAD2DEG)
    print('om ',nut.omega_moon % 360.0)
    # eps = 23.4387368 * constants.DEG2RAD
    eps = 23.4407685 * constants.DEG2RAD
    dpsi = -0.0034108 * constants.DEG2RAD
    omega_moon = 42.6046140 * constants.DEG2RAD
    Eq = rotations.equinox1982(dpsi, eps, omega_moon)
    gast = rotations.theta_GAST1982(gmst, Eq)
    assert_almost_equal(gast * constants.RAD2DEG, 312.8067654)


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



def test_appendix_c_conversion_from_TEME_to_ITRF_with_teme2ecef():
    """Test TEME to ITRF conversion

    References:
        Vallado et al., Revision 2
        Rhodes, Skyfield library, test_earth_satellites.py
    """
    rTEME = np.array([[5094.18016210], [6127.64465950], [6380.34453270]])
    
    # Apr 6, 2004,  07:51:28.386 UTC
    seconds = 28.386
    deltaUTC1 = -0.439961 # seconds
    s, us = np.divmod(seconds + deltaUTC1, 1)
    us *= 1e6  # microseconds
    dt = datetime.datetime(2004, 4, 6, 7, 51, int(s), int(us))
    jd = timefn.julian_date(dt)

    # Polar motion
    xp = -0.140682  # arcseconds
    yp = 0.333309  # arcseconds
    xp *= constants.ASEC2RAD
    yp *= constants.ASEC2RAD
    rITRF = rotations.teme2ecef(rTEME, jd, xp, yp)

    print(rITRF)
    assert_almost_equal(rITRF[0], -1033.47938300, decimal=4)
    assert_almost_equal(rITRF[1], 7901.29527540, decimal=4)
    assert_almost_equal(rITRF[2], 6380.35659580, decimal=4)


if __name__ == "__main__":
    import pytest
    pytest.main(['-v', __file__])