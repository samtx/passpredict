import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

from passpredict.constants import DEG2RAD, RAD2DEG, ASEC2RAD
from passpredict.utils import epoch_from_tle_datetime
from passpredict.rotations import transform
from passpredict.timefn import julian_date, jday2datetime


def test_eci2tod():
    """
    Appendix C, AIAA 2006-6753
    """
    rECI = np.array((3961.7442603, 6010.2156109, 4619.3625758))  # J2000, IAU 76-FK5
    xp = 0.098700 * ASEC2RAD
    yp = 0.286000 * ASEC2RAD
    dUT1 = 0.162360 # seconds
    dAT = 21  # seconds
    jd = julian_date(2000, 6, 28, 15, 8, 51.655000)
    jd_ut1 = 2451724.13115529340
    tt = 0.004904360547
    rTOD = np.array((3961.4214985, 6010.4782688, 4619.3015310))
    rPEF = np.array((298.8036328, -7192.3146229, 4619.3015310))
    rTEME = np.array((3961.0035498, 6010.7511740, 4619.3009301))


@pytest.mark.skip()
def test_ecef2eci():
    """
    Vallado, Eg. 3-15, p. 230
    """
    rECEF = np.array([
        [-1033.4793830],
        [ 7901.2952754],
        [ 6380.3565958]
    ])
    jdt_tt = 2453101.828154745
    jdt_ut1 = 2453101.827406783
    rECI = transform.ecef2eci(rECEF, jdt_ut1)
    print(f'dt_ut1     = {jday2datetime(jdt_ut1)}')
    print(f'jdt_ut1    = {jdt_ut1}')
    print(f'dt_tt      = {jday2datetime(jdt_tt)}')
    print(f'jdt_tt     = {jdt_tt}')
    assert_allclose(
        rECI,
        np.array([[5102.5096], [6123.01152], [6378.1363]]),
        atol=0.1  # need to make more accurate
    )


@pytest.mark.skip()
@pytest.mark.parametrize(
    'jd, rTEME, rECI_expected, xp, yp, decimal',
    [
        pytest.param(
            2453101.8274067827,   # Apr 6, 2004, 07:51:28.386 UTC,  -0.439961 sec deltaT
            np.array([[5094.18016210], [6127.64465950], [6380.34453270]]),
            np.array([-1033.47938300, 7901.29527540, 6380.35659580]),
            -0.140682 * ASEC2RAD,
            0.333309 * ASEC2RAD,
            4,
            id='Vallado, appendix c teme to itrf'
        ),
    ],
)
def test_teme2ecef_from_jd(jd, rTEME, rECI_expected, xp, yp, decimal):
    """
    Test TEME -> ITRF (ECEF) conversion
    """
    rECEF = transform.teme2ecef(rTEME, jd, xp, yp)
    assert_almost_equal(rECEF[0], rECI_expected[0], decimal=decimal)
    assert_almost_equal(rECEF[1], rECI_expected[1], decimal=decimal)
    assert_almost_equal(rECEF[2], rECI_expected[2], decimal=decimal)


@pytest.mark.skip()
@pytest.mark.parametrize(
    'jd, rTEME, rECI_expected, decimal',
    [
        pytest.param(
            2453101.8274067827,   # Apr 6, 2004, 07:51:28.386 UTC,  -0.439961 sec deltaT
            np.array((-9060.47373569, 4658.70952502, 813.68673153)),
            np.array((-9059.9413786, 4659.6972000, 813.9588875)),
            4,
            id='Vallado, appendix c teme to itrf'
        ),
    ],
)
def test_teme2eci_from_jd(jd, rTEME, rECI_expected, decimal):
    """
    Test TEME -> J2000 (ECI) conversion
    """
    rECI = transform.teme2eci(rTEME, jd)
    assert_almost_equal(rECI[0], rECI_expected[0], decimal=decimal)
    assert_almost_equal(rECI[1], rECI_expected[1], decimal=decimal)
    assert_almost_equal(rECI[2], rECI_expected[2], decimal=decimal)


@pytest.mark.skip()
@pytest.mark.parametrize(
    'epoch_string, rTEME, rECI_expected, decimal',
    [
        pytest.param(
            '00182.78495062',
            np.array((-9060.47373569, 4658.70952502, 813.68673153)),
            np.array((-9059.9413786, 4659.6972000, 813.9588875)),
            4, 
            id='Vallado, p. 234, epoch=00182.78495062'), # jd = 2451726.28495062
        pytest.param(
            '00179.78495062', 
            np.array((-9060.47373569, 4658.70952502, 813.68673153)),
            np.array((-9059.9510799, 4659.6807556, 813.9450451)),
            4, 
            id='Vallado, p. 234, epoch=00179.78495062'), # jd = 2451723.28495062
    ],
)
def test_teme2eci_from_epoch(epoch_string, rTEME, rECI_expected, decimal):
    """
    Test TEME -> J2000 (ECI) conversion

    Using 'of date'

    Reference:
        Vallado, p.234, Eq. 3-91
    """
    epoch = epoch_from_tle_datetime(epoch_string)
    jd = julian_date(epoch)
    rECI = transform.teme2eci(rTEME, jd)
    assert_almost_equal(rECI[0], rECI_expected[0], decimal=decimal)
    assert_almost_equal(rECI[1], rECI_expected[1], decimal=decimal)
    assert_almost_equal(rECI[2], rECI_expected[2], decimal=decimal)


if __name__ == "__main__":
    pytest.main(['-v', __file__])