import pytest
from numpy.testing import assert_almost_equal, assert_allclose
import numpy as np

import passpredict.rotations.nutation as nutation
from passpredict.constants import DEG2RAD, J2000, DJC, RAD2DEG
from passpredict.utils import epoch_from_tle_datetime
from passpredict.timefn import julian_date, jd2jc


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
    
    assert_almost_equal(nut80.M_moon*RAD2DEG, 314.9118590, decimal=5)
    assert_almost_equal(nut80.M_sun*RAD2DEG, 91.9379931, decimal=5)
    assert_almost_equal(nut80.u_M_moon*RAD2DEG, 169.0968272, decimal=5)
    assert_almost_equal(nut80.D_sun*RAD2DEG, 196.7518116, decimal=5)
    assert_almost_equal(nut80.omega_moon*RAD2DEG, 42.6046140, decimal=5)


@pytest.mark.skip()
@pytest.mark.parametrize(
    'tt, dpsi_deg, deps_deg, meaneps_deg',
    [
        pytest.param(0.004881175923888545, -0.004337544, -0.001247061, 23.43922657, id='Vallado, p.234, epoch=179.78495062'),
        pytest.param(0.004963311447502509, -0.004250260, -0.001260854, 23.43922657, id='Vallado, p.234, epoch=182.78495062'),
        pytest.param(        0.0426236319,   -0.0034108,    0.0020316, 23.4387368 , id='Vallado, p.230, eg.3-15')
    ]
)
def test_nut80_angles(tt, dpsi_deg, deps_deg, meaneps_deg):
    nut80 = nutation.nut80_fundamental_arguments(tt)
    dpsi, deps = nutation.nut80_angles(tt, nut80)
    meaneps = nutation.nut80_mean_eps(tt)
    assert_almost_equal(np.degrees(meaneps), meaneps_deg)
    assert_almost_equal(dpsi*RAD2DEG, dpsi_deg)
    assert_almost_equal(deps*RAD2DEG, deps_deg)


def test_nut80_angles_SOFA():
    """
    from SOFA validation routines
    http://www.iausofa.org/2019_0722_C/sofa/t_sofa_c.c
    t_nut80()
    """
    tt = 0.06
    nut80 = nutation.nut80_fundamental_arguments(tt)
    dpsi, deps = nutation.nut80_angles(tt, nut80)
    assert_almost_equal(dpsi, -0.9643658353226563966e-5, decimal=13)
    assert_almost_equal(deps,  0.4060051006879713322e-4, decimal=13)


@pytest.mark.skip()
@pytest.mark.parametrize(
    'tt, N_expected',
    [
        pytest.param(
            0.004881175923888545,
            np.array(
                (
                    ( 0.99999999955, 0.00000000000, 0.00003011190),
                    (-0.00000000066, 0.99999999976, 0.00002176637),
                    (-0.00003011190, 0.00002176637, 0.99999999931)
                )
            ), 
            id='Vallado, Appendix C, epoch=179.78495062'),
        pytest.param(
            0.004963311447502509,
            np.array(
                (
                    ( 0.99999999956, 0.00000000000, 0.00002950595),
                    (-0.00000000065, 0.99999999976, 0.00002200705),
                    (-0.00002950595, 0.00002200705, 0.99999999932)
                )
            ), 
            id='Vallado, Appendix C, epoch=182.78495062'),
    ]
)
def test_nutation_matrix(tt, N_expected):
    """
    Vallado, Appendix C
    """
    nut = nutation.fk5_nutation(tt)
    assert_allclose(nut.mtx, N_expected)



# jd1, jd2 = 2400000.5, 54388.0  # iauObl80
# meaneps = 0.4090751347643816218
# decimal = 14


if __name__ == "__main__":
    pytest.main(['-v', '--show-capture=all', __file__])

