# Test using pytest

import passpredictor as pp
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np

def test_site_declination_and_K():
    """
    Vallado, Eg 3-1
    """
    # Mt. Evans, Colorado
    phi_gd = 39.586667  # [deg]
    lmda = -105.640     # [deg]
    H_MSL = 4347.667    # [m]
    rdelta_calc, rK_calc = pp.site_declination_and_K(phi_gd, H_MSL)
    rdelta_true, rK_true = 4925.4298026, 4045.4937426
    assert_almost_equal(rdelta_calc, rdelta_true, decimal=2)
    assert_almost_equal(rK_calc, rK_true, decimal=2)

def test_site_declination_and_K_2():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007     # [deg]
    lmda = -104.883     # [deg]
    alt = 2187.    # [m]
    rdelta_calc, rK_calc = pp.site_declination_and_K(phi_gd, alt)
    rdelta_true, rK_true = 4964.5377, 3994.2955
    assert_almost_equal(rdelta_calc, rdelta_true, decimal=2)
    assert_almost_equal(rK_calc, rK_true, decimal=2)

def test_site_ECEF():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007     # [deg]
    lmda = -104.883     # [deg]
    alt = 2187.         # [m]
    r_site_ECEF_calc = pp.site_ECEF(phi_gd, lmda, alt)
    r_site_ECEF_true = np.array([-1275.1219, -4797.9890, 3994.2975], dtype=np.float64)
    assert_allclose(r_site_ECEF_calc, r_site_ECEF_true)


if __name__ == "__main__":
    test_site_declination_and_K_2()