# Test using pytest
from passpredict import kepler
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np
from datetime import datetime, timedelta


def test_coe2rv():
    """
    Vallado, Eg.2-6, p.116
    """
    p = 11067.790  # [km]
    e = 0.83285
    i = 87.87  # [deg]
    Omega = 227.89  # [deg]
    w = 53.38  # [deg]
    nu = 92.335  # [deg]
    r_calc, v_calc = kepler.coe2rv(p, e, i, Omega, w, nu)
    r_true = np.array([6525.344, 6861.535, 6449.125])  # [km]
    v_true = np.array([4.902276, 5.533124, -1.975709])  # [km/s]
    assert_allclose(r_calc, r_true, atol=1e-1)
    assert_allclose(v_calc, v_true, rtol=1e-5)


def test_rv2coe():
    """
    Vallado, Eg.2-5, p.114
    """
    rIJK = np.array([6524.834, 6862.875, 6448.296])  # [km]
    vIJK = np.array([4.901327, 5.533756, -1.976341])  # [km/s]
    p_true, a_true, e_true = 11067.790, 36127.343, 0.832853
    i_true, Omega_true, w_true = 87.870, 227.898, 53.38
    nu_true, u_true, lmda_true_true, what_true_true = (
        92.335,
        145.60549,
        55.282587,
        247.806,
    )
    coe = kepler.rv2coe(rIJK, vIJK, findall=True)
    assert_almost_equal(coe["p"], p_true, decimal=2)
    assert_almost_equal(coe["a"], a_true, decimal=2)
    assert_almost_equal(coe["e"], e_true, decimal=6)
    assert_almost_equal(coe["i"], i_true, decimal=3)
    assert_almost_equal(coe["Omega"], Omega_true, decimal=3)
    assert_almost_equal(coe["w"], w_true, decimal=2)
    assert_almost_equal(coe["nu"], nu_true, decimal=3)
    assert_almost_equal(coe["u"], u_true, decimal=0)
    assert_almost_equal(coe["lmda_true"], lmda_true_true, decimal=3)
    assert_almost_equal(coe["what_true"], what_true_true, decimal=3)


def test_rv2coe_to_coe2rv():
    rIJK = np.array([6524.834, 6862.875, 6448.296])  # [km]
    vIJK = np.array([4.901327, 5.533756, -1.976341])  # [km/s]
    coe = kepler.rv2coe(rIJK, vIJK, findall=True)
    p = coe.get("p")
    e = coe.get("e")
    i = coe.get("i")
    Omega = coe.get("Omega")
    w = coe.get("w")
    u = coe.get("u")
    lmda_true = coe.get("lmda_true")
    nu = coe.get("nu")
    what_true = coe.get("what_true")
    r_calc, v_calc = kepler.coe2rv(p, e, i, Omega, w, nu, u, lmda_true, what_true)
    assert_allclose(r_calc, rIJK)
    assert_allclose(v_calc, vIJK)


def test_coe2rv_to_rv2coe():
    p = 11067.790  # [km]
    e = 0.83285
    i = 87.87  # [deg]
    Omega = 227.89  # [deg]
    w = 53.38  # [deg]
    nu = 92.335  # [deg]
    r, v = kepler.coe2rv(p, e, i, Omega, w, nu)
    coeT = kepler.rv2coe(r, v, findall=True)
    assert_almost_equal(p, coeT["p"], decimal=2)
    assert_almost_equal(e, coeT["e"], decimal=6)
    assert_almost_equal(i, coeT["i"], decimal=3)
    assert_almost_equal(Omega, coeT["Omega"], decimal=3)
    assert_almost_equal(w, coeT["w"], decimal=2)
    assert_almost_equal(nu, coeT["nu"], decimal=3)


def test_nu2anomaly():
    """
    Vallado, p.87
    """
    nu = 60  # [deg]
    e = 0.999
    E_calc = kepler.nu2anomaly(nu, e)
    E_true = 1.4796584  # [deg]
    assert_almost_equal(E_calc, E_true)


def test_anomaly2nu():
    """
    Vallado, p.87
    """
    E = 1.4796584  # [deg]
    e = 0.999
    nu_calc = kepler.anomaly2nu(e, E)
    nu_true = 60.0  # [deg]
    assert_almost_equal(nu_calc, nu_true, decimal=6)


def test_kepEqtnE():
    """
    Vallado, Eg 2-1, p.66
    """
    M = 235.4  # [deg]
    e = 0.4
    E_calc = kepler.kepEqtnE(M, e)
    E_true = 220.512074767522  # [deg]
    assert_almost_equal(E_calc, E_true, decimal=12)


# def test_coe2rv_2():
#     """
#     Vallado, Eg.11-6, COE from 2.4.2
#     """
#     a = 6768.3568   # [km]
#     e = 0.0005770
#     p = a*(1-e**2)
#     i = 51.6190       # [deg]
#     Omega = 13.3340  # [deg]
#     w = 102.5680       # [deg]
#     M = 257.5950  # mean anomaly
#     E = kepler.kepEqtnE(M, e)
#     nu = kepler.anomaly2nu(e, E)
#     r_calc, v_calc = kepler.coe2rv(p, e, i, Omega, w, nu)
#     r_true = np.array([6585.038266, 1568.184321, 9.116355])  # [km]
#     v_true = np.array([-1.1157766, 4.6316816, 6.0149576]) # [km/s]
#     assert_allclose(r_calc, r_true, atol=1e-1, verbose=True)
#     assert_allclose(v_calc, v_true, rtol=1e-5, verbose=True)


# def test_pkepler():
#     """
#     Vallado, Eg11-6
#     """
#     r0 = np.array([6585.038266, 1568.184321, 9.116355])
#     v0 = np.array([-1.1157766, 4.6316816, 6.0149576])
#     jd0 = 2450540.400
#     ndt = 7.889e-5
#     # jdf = kepler.julian_date(1997, 4, 2, 13, 8, 0)
#     jdf = 2450540.5472
#     dt = jdf - jd0
#     rf, vf = kepler.pkepler_rv(r0, v0, dt, ndt)
#     rf_true = np.array([-2811.2769, 3486.2632, 5069.5763])
#     vf_true = np.array([-6.859691, -2.964792, -1.764721])
#     assert_allclose(rf, rf_true)
#     assert_allclose(vf, vf_true)


# def test_pkepler_2():
#     """
#     From Fortran test routines
#     """

#     r0 = np.array([1.02, 1.03, 1.04])
#     v0 = np.array([0.784, 0.27, 0.012])
#     ndt, nddt, dtsec = 0.000003, 0.00000005, 23.43
