import pytest
from numpy.testing import assert_almost_equal

import passpredict.nutation as nutation
from passpredict.constants import DEG2RAD, J2000, DJC, RAD2DEG

def test_nut80_fundamental_arguments():
    """
    Vallado, Eg. 3-15, p.230
    """
    tt = 0.0426236319  # April 6, 2004, 07:51:28.386009, dUT1=-0.4399619, dAT=32
    nut80 = nutation.nut80_fundamental_arguments(tt)
    print(nut80)
    # assert_almost_equal(nut80.M_moon*RAD2DEG, 314.9118590)
    # assert_almost_equal(nut80.M_sun*RAD2DEG, 91.9379931)
    # assert_almost_equal(nut80.u_M_moon*RAD2DEG, 169.0968272)
    # assert_almost_equal(nut80.D_sun*RAD2DEG, 196.7518116)
    # assert_almost_equal(nut80.omega_moon*RAD2DEG, 42.6046140)
    
    assert_almost_equal(nut80.M_moon, 314.9118590, decimal=5)
    assert_almost_equal(nut80.M_sun, 91.9379931, decimal=5)
    assert_almost_equal(nut80.u_M_moon, 169.0968272, decimal=5)
    assert_almost_equal(nut80.D_sun, 196.7518116, decimal=5)
    assert_almost_equal(nut80.omega_moon, 42.6046140, decimal=5)


def test_nut80_angles():
    """
    Vallado, Eg. 3-15, p.230
    """
    tt = 0.0426236319  # April 6, 2004, 07:51:28.386009, dUT1=-0.4399619, dAT=32
    nut80 = nutation.nut80_fundamental_arguments(tt)
    dpsi, deps = nutation.nut80_angles(tt, nut80)
    print(f'dpsi={dpsi:.8f},  deps={deps:.8f}')
    assert_almost_equal(dpsi*RAD2DEG, -0.0034108)
    assert_almost_equal(deps*RAD2DEG, 0.0020316)


def test_nut80_angles_SOFA():
    """
    from SOFA validation routines
    http://www.iausofa.org/2019_0722_C/sofa/t_sofa_c.c
    t_nut80()
    """

    jd1 = 2400000.5
    jd2 = 53736.0
    tt = ((jd1 - J2000) + jd2)/ DJC
    print(f'tt = {tt:.12f}')
    nut80 = nutation.nut80_fundamental_arguments(tt)
    dpsi, deps = nutation.nut80_angles(tt, nut80)
    print(f'dpsi={dpsi:.8f},  deps={deps:.8f}')
    assert_almost_equal(dpsi, -0.9643658353226563966e-5, decimal=13)
    assert_almost_equal(deps,  0.4060051006879713322e-4, decimal=13)


def test_nut80_mean_eps():
    """
    Vallado, Eg. 3-15, p.230
    """
    tt = 0.0426236319  # April 6, 2004, 07:51:28.386009, dUT1=-0.4399619, dAT=32
    meaneps = nutation.nut80_mean_eps(tt)
    meaneps *= RAD2DEG
    print(f'meaneps = {meaneps:0.8f}')
    assert_almost_equal(meaneps, 23.4387368)



if __name__ == "__main__":
    pytest.main(['-v', '--show-capture=all', __file__])

