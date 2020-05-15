# Test using pytest

from passpredictor import predict
from passpredictor import timefn
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np
import datetime
from passpredictor.constants import ASEC2RAD

def test_satellite_visible():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([[-2811.2769, 3486.2632, 5069.5763]]).T  # ECI coords
    rsite = np.array([[-3414.0283, 3258.1636, 4276.1212]]).T  # ECI coords
    rho = np.array([[-773.8654, -581.4980, 328.8145]]).T  # SEZ coords
    dt = datetime.datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = np.array([timefn.julian_date2(dt)])
    vis = predict.satellite_visible(rsat, rsite, rho, jdt)
    assert vis[0] > 2


# def test_riseset():
#     """Eg. ,
#     Test orbit propagation
#     """
#         # Obj, n [rev/solar day],         e,  i [deg],  w, Omega, M
#     orbital_objects = [
#         (   1,        1.00272141, 0.0000032,  00.0956, 0.,    0., 0.),
#         (   2,        8.36589235, 0.0080158,  90.0175, 0.,    0., 0.),
#         (   3,        0.24891961, 0.9363060,  64.9874, 0.,    0., 0.),
#         (   4,        0.21467209, 0.0668128,  57.3500, 0.,    0., 0.),
#         (   5,       13.37659679, 0.0145072,  90.2619, 0.,    0., 0.),
#         (   6,       16.09769232, 0.0078742,  82.8709, 0.,    0., 0.),
#         (   7,        1.00271920, 0.0003109,  00.0099, 0.,    0., 0.),
#         (   8,       12.41552416, 0.0036498,  74.0186, 0.,    0., 0.),
#         (   9,       13.84150848, 0.0048964, 144.6414, 0.,    0., 0.)
#     ]


def test_apredictendix_c_conversion_from_TEME_to_ITRF_UTC1():
    """Test TEME to ITRF conversion

    References:
        Vallado et al., Revision 2
        Rhodes, Skyfield library, test_earth_satellites.py
    """
    rTEME = np.array([5094.18016210, 6127.64465950, 6380.34453270])
    vTEME = np.array([-4.746131487, 0.785818041, 5.531931288])
    vTEME = vTEME * 24.0 * 60.0 * 60.0  # km/s to km/day

    # Apr 6, 2004,  07:51:28.386 UTC
    # deltaUTC1 = -0.439961 seconds
    ms = int(1e6) + 386000 - 439961
    dt = datetime.datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)

    # Polar motion
    xp = -0.140682  # arcseconds
    yp = 0.333309  # arcseconds
    xp *= ASEC2RAD
    yp *= ASEC2RAD
    # xp = yp = 0.
    rITRF, vITRF = predict.TEME_to_ITRF(jd, rTEME, vTEME, xp, yp)

    print(rITRF)
    assert_almost_equal(rITRF[0], -1033.47938300, decimal=4)
    assert_almost_equal(rITRF[1], 7901.29527540, decimal=4)
    assert_almost_equal(rITRF[2], 6380.35659580, decimal=4)
