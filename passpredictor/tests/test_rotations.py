# test rotations.py
from .. import rotations
from .. import predict
from .. import constants
from .. import timefn
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np
from datetime import datetime, timedelta


def test_ECEF_to_SEZ():
    """
    Vallado, Eg. 11-6, p.912
    """
    phi = 42.38    # latitude, deg
    lmda = -71.13  # longitude, deg
    # lmda = 136.2944
    h = 24         # height, m
    rsat = np.array([885.7296, -4389.3856, 5070.1765])
    rsite = predict.site_ECEF2(phi, lmda, h)
    rhoECEF = rsat - rsite
    print(rhoECEF)
    rSEZ = rotations.ECEF_to_SEZ(rhoECEF, phi, lmda)
    rSEZ_true = np.array([-773.8654, -581.4980, 328.8145])
    np.set_printoptions(precision=8)
    # print(rSEZ)
    for i in [0, 1, 2]:
        assert_almost_equal(rSEZ[i], rSEZ_true[i], decimal=0, verbose=True)


def test_fk5_precession():
    """
    Vallado, Eg. 3-15, p.231
    """
    # April 6, 2004, 07:51:28.386009 UTC
    tt = 0.0426236319  # Julian centuries since J2000
    zeta, theta, z = rotations.fk5_precession(tt)
    assert_almost_equal(zeta, 0.0273055*constants.DEG2RAD, decimal=9)
    assert_almost_equal(theta, 0.0237306*constants.DEG2RAD, decimal=9)
    assert_almost_equal(z, 0.0273059*constants.DEG2RAD, decimal=9)


def test_precess_rotation():
    """
    Vallado, Eg. 3-15, p.231
    """
    rMOD = np.array([5094.0283745, 6127.8708164, 6380.2485164])
    zeta = 0.0273055 * predict.DEG2RAD
    theta = 0.0237306 * predict.DEG2RAD
    z = 0.0273059 * predict.DEG2RAD
    rGCRF = rotations.precess_rotation(rMOD, zeta, theta, z)
    rGCRF_true = np.array([5102.508958, 6123.011401, 6378.136928])
    assert_allclose(rGCRF, rGCRF_true)


def test_theta_GMST1982():
    """Compute the Greenwich Mean Sidereal Time

    References:
        Vallado, Eg 3-5, p.188
    """
    dt = datetime(1992, 8, 20, 12, 14, 0)  # Aug 20, 1992, 12:14 PM UT1
    jd = timefn.julian_date(dt)
    theta, thetadt = rotations.theta_GMST1982(jd)
    theta *= predict.RAD2DEG  # convert from radians to degrees
    assert_almost_equal(theta, 152.578787810, decimal=9)  # degrees
    # assert_almost_equal(thetadt, 152.578787810)  # degrees


def test_theta_GMST1982_2():
    """Compute the Greenwich Mean Sidereal Time

    References
        Vallado, Eg. 3-15, p.230
    """
    # April 6, 2004, 07:51:28.386 009UTC, not UTC1
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)
    theta, thetadt = rotations.theta_GMST1982(jd)
    theta *= predict.RAD2DEG
    assert_almost_equal(theta, 312.8098943, decimal=6)


def test_IJK2SEZ():
    """
    Curtis, Eg. 5.9, p.270
    """
    from math import sin, cos, radians, asin, acos, degrees
    # Satellite position
    r = np.array([-2032.4, 4591.2, -4544.8])  # km, geocentric equatorial pos.
    # Location position
    H = 0   # elevation, sea level
    phi = -40    # deg latitude
    theta = 110  # deg, local sidereal time
    R_obs = np.array([-1673.0, 4598.0, -4078.0])  # km, from Eq. 5.56

    # sez = predict.IJK2SEZ(r, phi, theta, H)
    rho = r - R_obs  # relative position vector
    phi_rad = radians(phi)
    theta_rad = radians(theta)

    # rotation matrix, from Eq. 5.62a
    Q = np.array([
       [-sin(theta_rad), cos(theta_rad), 0.],
       [-sin(phi_rad)*cos(theta_rad), -sin(phi_rad)*sin(theta_rad), cos(phi_rad)],
       [cos(phi_rad)*cos(theta_rad), cos(phi_rad)*sin(theta_rad), sin(phi_rad)]
    ])

    rho2 = np.dot(Q, rho)

    rho2_hat = (1/np.linalg.norm(rho))*rho2

    a = degrees(asin(rho2_hat[2]))
    A = degrees(acos(rho2_hat[1]/cos(radians(a))))

    A_true = 129.8  # deg, azimuth
    a_true = 41.41  # deg, angular elevation
    assert_almost_equal(A, A_true, decimal=1)
    assert_almost_equal(a, a_true, decimal=1)


def test_thetaGMST_from_skyfield():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)
    theta, thetadt = rotations.theta_GMST1982(jd)
    assert_almost_equal(theta, 5.459562584754709, decimal=15)
