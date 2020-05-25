# Test using pytest
from passpredict import topocentric
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np


def test_site_declination_and_K():
    """
    Vallado, Eg 3-1
    """
    # Mt. Evans, Colorado
    phi_gd = 39.586667  # [deg]
    H_MSL = 4347.667  # [m]
    rdelta_calc, rK_calc = topocentric.site_declination_and_K(phi_gd, H_MSL)
    rdelta_true, rK_true = 4925.4298026, 4045.4937426
    assert_almost_equal(rdelta_calc, rdelta_true, decimal=3)
    assert_almost_equal(rK_calc, rK_true, decimal=3)


def test_site_declination_and_K_2():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007  # [deg]
    alt = 2187.0  # [m]
    rdelta_calc, rK_calc = topocentric.site_declination_and_K(phi_gd, alt)
    rdelta_true, rK_true = 4964.5377, 3994.2955
    assert_almost_equal(rdelta_calc, rdelta_true, decimal=2)
    assert_almost_equal(rK_calc, rK_true, decimal=2)


def test_sgp4_fk5():
    """
    Vallado, p. 233-234
    """
    pass


def test_site_ECEF():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007  # [deg]
    lmda = -104.883  # [deg]
    alt = 2187.0  # [m]
    r_ECEF = topocentric.site_ECEF(phi_gd, lmda, alt)
    r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    for i in [0, 1, 2]:
        assert_almost_equal(r_ECEF[i], r_ECEFtrue[i], decimal=4)


def test_site_ECEF2():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007  # [deg]
    lmda = -104.883  # [deg]
    alt = 2187.0  # [m]
    r_ECEF = topocentric.site_ECEF2(phi_gd, lmda, alt)
    r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    for i in [0, 1, 2]:
        assert_almost_equal(r_ECEF[i], r_ECEFtrue[i], decimal=4)


def test_site_ECEF2_v2():
    phi = 42.38  # latitude, deg
    lmda = -71.13  # longitude, deg
    h = 24  # height, m
    rsite = topocentric.site_ECEF2(phi, lmda, h)
    rtrue = np.array([1526.122, -4465.064, 4276.894])
    print(rsite)
    for i in [0, 1, 2]:
        assert_almost_equal(rsite[i], rtrue[i], decimal=3, verbose=True)
