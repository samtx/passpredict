# Test using pytest

from .. import predict
from .. import timefn
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np
from datetime import datetime, timedelta
from ..constants import (
    R_EARTH,
    R2_EARTH,
    e_EARTH,
    e2_EARTH,
    MU,
    J2,
    J2000,
    AU_M,
    AU_KM,
    ASEC360,
    DAY_S,
    ASEC2RAD,
    DEG2RAD,
    RAD2DEG,
    tau,
)


def test_site_declination_and_K():
    """
    Vallado, Eg 3-1
    """
    # Mt. Evans, Colorado
    phi_gd = 39.586667  # [deg]
    H_MSL = 4347.667  # [m]
    rdelta_calc, rK_calc = predict.site_declination_and_K(phi_gd, H_MSL)
    rdelta_true, rK_true = 4925.4298026, 4045.4937426
    assert_almost_equal(rdelta_calc, rdelta_true, decimal=3)
    assert_almost_equal(rK_calc, rK_true, decimal=3)


def test_site_declination_and_K_2():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007  # [deg]
    alt = 2187.0  # [m]
    rdelta_calc, rK_calc = predict.site_declination_and_K(phi_gd, alt)
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
    r_ECEF = predict.site_ECEF(phi_gd, lmda, alt)
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
    r_ECEF = predict.site_ECEF2(phi_gd, lmda, alt)
    r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    for i in [0, 1, 2]:
        assert_almost_equal(r_ECEF[i], r_ECEFtrue[i], decimal=4)


def test_site_ECEF2_v2():
    phi = 42.38  # latitude, deg
    lmda = -71.13  # longitude, deg
    h = 24  # height, m
    rsite = predict.site_ECEF2(phi, lmda, h)
    rtrue = np.array([1526.122, -4465.064, 4276.894])
    print(rsite)
    for i in [0, 1, 2]:
        assert_almost_equal(rsite[i], rtrue[i], decimal=3, verbose=True)


def test_sun_pos():
    """
    Vallado, Eg. 5-1, p. 280
    """
    dt = datetime(2006, 4, 2)  # April 2, 2006, 00:00 UTC
    jdt = timefn.julian_date2(dt)
    assert_almost_equal(jdt, 2453827.5, decimal=12)
    jdt = np.asarray(jdt)
    r = predict.sun_pos(jdt)
    r_true = np.array([146186212, 28788976, 12481064], dtype=np.float)
    r_true = np.reshape(r_true, (3, 1))
    assert_allclose(r, r_true, rtol=1e-4)


def test_sun_pos_2():
    """
    Vallado, Eg. 11-6, p. 913
    """
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = timefn.julian_date2(dt)
    jdt = np.asarray(jdt)
    r = predict.sun_pos(jdt) / AU_KM
    r_true = np.array([0.9765, 0.1960, 0.0850], dtype=np.float)
    r_true = np.reshape(r_true, (3, 1))
    assert_allclose(r, r_true, rtol=1e-3)


def test_sun_sat_angle():
    """
    Vallado, Eg. 11-6, p.913
    """
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = timefn.julian_date2(dt)
    jdt = np.asarray(jdt)
    rsun = predict.sun_pos(jdt)
    rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
    rsun = np.atleast_2d(rsun)
    rsat = np.atleast_2d(rsat).T
    sunangle = predict.sun_sat_angle(rsat, rsun) * RAD2DEG
    assert_almost_equal(sunangle, 76.0407, decimal=3)


def test_sun_sat_angle2():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
    rsun = np.array([0.9765, 0.1960, 0.0850]) * AU_KM
    sunangle = predict.sun_sat_angle(rsat, rsun) * RAD2DEG
    assert_almost_equal(sunangle, 76.0407, decimal=3)


def test_satellite_visible():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([[-2811.2769, 3486.2632, 5069.5763]]).T  # ECI coords
    rsite = np.array([[-3414.0283, 3258.1636, 4276.1212]]).T  # ECI coords
    rho = np.array([[-773.8654, -581.4980, 328.8145]]).T  # SEZ coords
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
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


def test_sun_sat_orthogonal_distance():
    """
    Vallado, Eg. 11-6, p.913
    """
    r = np.array([-2811.2769, 3486.2632, 5069.5763])  # sat, ECI coordinates
    zeta = 76.0407  # deg
    dist = predict.sun_sat_orthogonal_distance(r, zeta * predict.DEG2RAD)
    assert_almost_equal(dist, 6564.6870, decimal=4)


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
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
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


# def test_sgp4():
#     """
#     Vallado, Eg.11-6 using SGP4 as propagator
#     """
#     tle1 = '1 16609U 86017A   93352.53502934  .00007889  0000  0  10529-3 0   342'
#     tle2 = '2 16609  51.6190  13.3340 0005770 102.5680 257.5950 15.59114070447869'
#     sat = predict.twoline2rv(tle1, tle2, predict.wgs72)
#     from predictrint import predictrint
#     predictrint(sat)
#     r0T = np.array([6585.038266, 1568.184321, 9.116355])
#     v0T = np.array([-1.1157766, 4.6316816, 6.0149576])
#     jd0T = 2450540.400
#     assert_almost_equal(predict.julian_date(sat.epoch),jd0T)


# def test_IJK2SEZ_2():
#     """
#     Vallado, p. 164
#     """

#     # example from Vallado, p.913
#     phi = 42.38    # degrees latitude
#     lmda = -71.13  # degrees local sidereal time
#     h_ellp = 24

#     r_sat_ECEF = np.array([885.7296, -4389.3856, 5069.5763])  # km
#     r_site_ECEF = predict.site_ECEF(phi, lmda, h_ellp)

#     from math import radians, degrees, sin, cos
#     # Rotation matrix IJK -> SEZ
#     angle2 = radians(90-phi)
#     angle3 = radians(lmda)

#     R = np.matmul(predict.rot2(angle2), predict.rot3(angle3))

#     # Get position vector from site to satellite
#     rho_ECEF = r_sat_ECEF - r_site_ECEF

#     # Rotate position vector to SEZ coordinates
#     rho_SEZ = np.matmul(R, rho_ECEF)

#     rho_SEZ_exact = np.array([-773.8654, -581.4980, 328.8145])  # km

#     assert_allclose(rho_SEZ, rho_SEZ_exact)

if __name__ == "__main__":
    # test_julian_date_vectorized()
    # test_theta_GMST1982()
    # test_apredictendix_c_conversion_from_TEME_to_ITRF_UTC1()
    # test_jd_from_skyfield2()
    # test_jd_from_skyfield3()
    # test_site_declination_and_K()
    # test_site_ECEF2()
    # test_ECEF_to_SEZ()
    #     test_satellite_visible()

    def azm(s, e):
        out = np.arctan2(s, e) * 180 / np.pi + 90
        print(out)
        if s < 0 and e < 0:
            out = out % 360
        return out

    s, e = -1, 1
    print(azm(s, e))
