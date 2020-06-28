import pytest
from numpy.testing import assert_almost_equal

import passpredict.nutation as nutation
from passpredict.constants import DEG2RAD, RAD2DEG

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


if __name__ == "__main__":
    pytest.main(['-v', __file__])
