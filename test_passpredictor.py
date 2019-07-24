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
    assert_almost_equal(rdelta_calc, rdelta_true, decimal=3)
    assert_almost_equal(rK_calc, rK_true, decimal=3)

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


def test_coe2rv():
    """
    Vallado, Eg.2-6, p.116
    """
    p = 11067.790   # [km]
    e = 0.83285
    i = 87.87       # [deg]
    Omega = 227.89  # [deg]
    w = 53.38       # [deg]
    nu = 92.335     # [deg]
    r_calc, v_calc = pp.coe2rv(p, e, i, Omega, w, nu)
    r_true = np.array([6525.344, 6861.535, 6449.125])  # [km] 
    v_true = np.array([4.902276, 5.533124, -1.975709]) # [km/s]
    assert_allclose(r_calc, r_true, atol=1e-1)
    assert_allclose(v_calc, v_true, rtol=1e-5)


def test_pkepler():
    """ 
    Test orbit propagation
    """
        # Obj, n [rev/solar day],         e,  i [deg],  w, Omega, M
    orbital_objects = [
        (   1,        1.00272141, 0.0000032,  00.0956, 0.,    0., 0.),
        (   2,        8.36589235, 0.0080158,  90.0175, 0.,    0., 0.),
        (   3,        0.24891961, 0.9363060,  64.9874, 0.,    0., 0.),
        (   4,        0.21467209, 0.0668128,  57.3500, 0.,    0., 0.),
        (   5,       13.37659679, 0.0145072,  90.2619, 0.,    0., 0.),
        (   6,       16.09769232, 0.0078742,  82.8709, 0.,    0., 0.),
        (   7,        1.00271920, 0.0003109,  00.0099, 0.,    0., 0.),
        (   8,       12.41552416, 0.0036498,  74.0186, 0.,    0., 0.),
        (   9,       13.84150848, 0.0048964, 144.6414, 0.,    0., 0.)
    ] 

def test_nu2anomaly():
    """
    Vallado, p.87
    """
    nu = 60   # [deg]
    e = 0.999
    E_calc = pp.nu2anomaly(nu, e)
    E_true = 1.4796584  # [deg]
    assert_almost_equal(E_calc, E_true)


def test_anomaly2nu():
    """
    Vallado, p.87
    """
    E = 1.4796584  # [deg]
    e = 0.999
    nu_calc = pp.anomaly2nu(e, E)
    nu_true = 60.  # [deg]
    assert_almost_equal(nu_calc, nu_true, decimal=6)


def test_kepEqtnE():
    """
    Vallado, Eg 2-1, p.66
    """
    M = 235.4 # [deg]
    e = 0.4
    E_calc = pp.kepEqtnE(M, e)
    E_true = 220.512074767522 # [deg]
    assert_almost_equal(E_calc, E_true)

if __name__ == "__main__":
    test_coe2rv()