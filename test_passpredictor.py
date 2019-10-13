# Test using pytest

import passpredictor as pp
from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np
from datetime import datetime, timedelta


def test_site_declination_and_K():
    """
    Vallado, Eg 3-1
    """
    # Mt. Evans, Colorado
    phi_gd = 39.586667  # [deg]
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
    alt = 2187.    # [m]
    rdelta_calc, rK_calc = pp.site_declination_and_K(phi_gd, alt)
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
    phi_gd = 39.007     # [deg]
    lmda = -104.883     # [deg]
    alt = 2187.         # [m]
    r_ECEF = pp.site_ECEF(phi_gd, lmda, alt)
    r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    for i in [0, 1, 2]:
        assert_almost_equal(r_ECEF[i], r_ECEFtrue[i], decimal=4)


def test_site_ECEF2():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007     # [deg]
    lmda = -104.883     # [deg]
    alt = 2187.         # [m]
    r_ECEF = pp.site_ECEF2(phi_gd, lmda, alt)
    r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    for i in [0, 1, 2]:
        assert_almost_equal(r_ECEF[i], r_ECEFtrue[i], decimal=4)


def test_site_ECEF2_v2():
    phi = 42.38    # latitude, deg
    lmda = -71.13  # longitude, deg
    h = 24         # height, m
    rsite = pp.site_ECEF2(phi, lmda, h)
    rtrue = np.array([1526.122, -4465.064, 4276.894])
    print(rsite)
    for i in [0, 1, 2]:
        assert_almost_equal(rsite[i], rtrue[i], decimal=3, verbose=True)


def test_ECEF_to_SEZ():
    """
    Vallado, Eg. 11-6, p.912
    """
    phi = 42.38    # latitude, deg
    lmda = -71.13  # longitude, deg
    # lmda = 136.2944
    h = 24         # height, m
    rsat = np.array([885.7296, -4389.3856, 5070.1765])
    rsite = pp.site_ECEF2(phi, lmda, h)
    rhoECEF = rsat - rsite
    print(rhoECEF)
    rSEZ = pp.ECEF_to_SEZ(rhoECEF, phi, lmda)
    rSEZ_true = np.array([-773.8654, -581.4980, 328.8145])
    np.set_printoptions(precision=8)
    # print(rSEZ)
    for i in [0, 1, 2]:
        assert_almost_equal(rSEZ[i], rSEZ_true[i], decimal=0, verbose=True)


def test_sun_pos():
    """
    Vallado, Eg. 5-1, p. 280
    """
    dt = datetime(2006, 4, 2)  # April 2, 2006, 00:00 UTC
    jdt = pp.julian_date2(dt)
    assert_almost_equal(jdt, 2453827.5, decimal=12)
    jdt = np.asarray(jdt)
    r = pp.sun_pos(jdt)
    r_true = np.array([146186212, 28788976, 12481064], dtype=np.float)
    r_true = np.reshape(r_true, (3, 1))
    assert_allclose(r, r_true, rtol=1e-4)


def test_sun_pos_2():
    """
    Vallado, Eg. 11-6, p. 913
    """
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = pp.julian_date2(dt)
    jdt = np.asarray(jdt)
    r = pp.sun_pos(jdt) / pp.AU_KM
    r_true = np.array([0.9765, 0.1960, 0.0850], dtype=np.float)
    r_true = np.reshape(r_true, (3, 1))
    assert_allclose(r, r_true, rtol=1e-3)


def test_sun_sat_angle():
    """
    Vallado, Eg. 11-6, p.913
    """
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = pp.julian_date2(dt)
    jdt = np.asarray(jdt)
    rsun = pp.sun_pos(jdt)
    rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
    rsun = np.atleast_2d(rsun)
    rsat = np.atleast_2d(rsat).T
    sunangle = pp.sun_sat_angle(rsat, rsun) * pp.RAD2DEG
    assert_almost_equal(sunangle, 76.0407, decimal=3)


def test_sun_sat_angle2():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
    rsun = np.array([0.9765, 0.1960, 0.0850]) * pp.AU_KM
    sunangle = pp.sun_sat_angle(rsat, rsun) * pp.RAD2DEG
    assert_almost_equal(sunangle, 76.0407, decimal=3)


def test_satellite_visible():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([[-2811.2769, 3486.2632, 5069.5763]]).T  # ECI coords
    rsite = np.array([[-3414.0283, 3258.1636, 4276.1212]]).T      # ECI coords
    rho = np.array([[-773.8654, -581.4980, 328.8145]]).T   # SEZ coords
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = np.array([pp.julian_date2(dt)])
    vis = pp.satellite_visible(rsat, rsite, rho, jdt )
    assert(vis[0] > 2)


def test_fk5_precession():
    """
    Vallado, Eg. 3-15, p.231
    """
    # April 6, 2004, 07:51:28.386009 UTC
    tt = 0.0426236319  # Julian centuries since J2000
    zeta, theta, z = pp.fk5_precession(tt)
    assert_almost_equal(zeta, 0.0273055*pp.DEG2RAD, decimal=9)
    assert_almost_equal(theta, 0.0237306*pp.DEG2RAD, decimal=9)
    assert_almost_equal(z, 0.0273059*pp.DEG2RAD, decimal=9)


def test_precess_rotation():
    """
    Vallado, Eg. 3-15, p.231
    """
    rMOD = np.array([5094.0283745, 6127.8708164, 6380.2485164])
    zeta = 0.0273055 * pp.DEG2RAD
    theta = 0.0237306 * pp.DEG2RAD
    z = 0.0273059 * pp.DEG2RAD
    rGCRF = pp.precess_rotation(rMOD, zeta, theta, z)
    rGCRF_true = np.array([5102.508958, 6123.011401, 6378.136928])
    assert_allclose(rGCRF, rGCRF_true)


def test_sun_sat_orthogonal_distance():
    """
    Vallado, Eg. 11-6, p.913
    """
    r = np.array([-2811.2769, 3486.2632, 5069.5763])  # sat, ECI coordinates
    zeta = 76.0407  # deg
    dist = pp.sun_sat_orthogonal_distance(r, zeta * pp.DEG2RAD)
    assert_almost_equal(dist, 6564.6870, decimal=4)


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


def test_utc2tt():
    """Vallado, Eg. 3-7"""
    # Mountain standard time (UTC-6)
    dt = datetime(2004, 5, 14, 10, 43, 0)
    deltaAT = 32  # sec
    deltaUT1 = -0.463326  # sec
    # get UTC by adjusting from MST timezone
    dt += timedelta(hours=6)
    jd_utc = pp.julian_date(dt)
    tt = pp.utc2tt(jd_utc, deltaAT=deltaAT, deltaUT1=deltaUT1)
    assert_almost_equal(tt, 0.043674121031, decimal=12)


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
    M = 235.4  # [deg]
    e = 0.4
    E_calc = pp.kepEqtnE(M, e)
    E_true = 220.512074767522  # [deg]
    assert_almost_equal(E_calc, E_true, decimal=12)


def test_julian_date():
    """
    Vallado, eg.3-4
    """
    yr, mo, dy = 1996, 10., 26.
    hr, mn, sec = 14., 20., 0.
    jd = pp.julian_date(yr, mo, dy, hr, mn, sec)
    jdT = 2450383.09722222
    assert_almost_equal(jd, jdT, decimal=8)


def test_julian_date_datetime():
    """
    Vallado, eg.3-4
    """
    yr, mo, dy = 1996, 10, 26
    hr, mn, sec = 14, 20, 0
    dt = datetime(yr, mo, dy, hr, mn, sec)
    jd = pp.julian_date(dt)
    jdT = 2450383.09722222
    assert_almost_equal(jd, jdT, decimal=8)


# def test_julian_date_datetime2():
#     """
#     Vallado, eg. 3-15, p. 230
#     """
#     dt = datetime(2004, 4, 6, 7, 51, )
#     jd = pp.julian_date(dt)
#     jdT = 2450383.09722222
#     assert_almost_equal(jd, jdT, decimal=8)


def test_julian_date_vectorized():
    """Use an array of datetimes to find the Julian Date"""
    dt_ary = np.arange('2019-09-14T00:00:00', '2019-10-07T00:00:00', 200, dtype='datetime64')
    jd_vectorized = np.vectorize(pp.julian_date)
    jd_ary = jd_vectorized(dt_ary)
    # print(jd_ary)


def test_theta_GMST1982():
    """Compute the Greenwich Mean Sidereal Time

    References:
        Vallado, Eg 3-5, p.188
    """
    dt = datetime(1992, 8, 20, 12, 14, 0)  # Aug 20, 1992, 12:14 PM UT1
    jd = pp.julian_date(dt)
    theta, thetadt = pp.theta_GMST1982(jd)
    theta *= pp.RAD2DEG  # convert from radians to degrees
    assert_almost_equal(theta, 152.578787810, decimal=9)  # degrees
    # assert_almost_equal(thetadt, 152.578787810)  # degrees


def test_theta_GMST1982_2():
    """Compute the Greenwich Mean Sidereal Time

    References
        Vallado, Eg. 3-15, p.230
    """
    # April 6, 2004, 07:51:28.386 009UTC, not UTC1
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = pp.julian_date(dt)
    theta, thetadt = pp.theta_GMST1982(jd)
    theta *= pp.RAD2DEG
    assert_almost_equal(theta, 312.8098943, decimal=6)


def test_appendix_c_conversion_from_TEME_to_ITRF_UTC1():
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
    jd = pp.julian_date(dt)

    # Polar motion
    xp = -0.140682  # arcseconds
    yp = 0.333309   # arcseconds
    xp *= pp.ASEC2RAD
    yp *= pp.ASEC2RAD
    # xp = yp = 0.
    rITRF, vITRF = pp.TEME_to_ITRF(jd, rTEME, vTEME, xp, yp)

    print(rITRF)
    assert_almost_equal(rITRF[0], -1033.47938300, decimal=4)
    assert_almost_equal(rITRF[1], 7901.29527540, decimal=4)
    assert_almost_equal(rITRF[2], 6380.35659580, decimal=4)


def test_jd_from_skyfield():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = pp.julian_date(dt)
    assert_almost_equal(jd, 2453101.8274067827, decimal=12)


def test_jd_from_skyfield2():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = pp.julian_date2(dt)
    jd_desired = 2453101.8274067827
    print(f'jd   ={jd:20.15f}')
    print(f'jddes={jd_desired:20.15f}')
    print(f'diff ={jd-jd_desired:0.15f}')
    assert_almost_equal(jd, 2453101.8274067827, decimal=12)


def test_jd_from_skyfield3():
    """From skyfield.tests.test_earth_satellites.py"""
    sec = 28.386 - 0.439961
    yr, mo, dy = 2004, 4, 6
    hr, mn = 7, 51
    jd = pp.julian_date2(yr, mo, dy, hr, mn, sec)
    jd_desired = 2453101.8274067827
    print(f'jd   ={jd:20.15f}')
    print(f'jddes={jd_desired:20.15f}')
    print(f'diff ={jd-jd_desired:0.15f}')
    assert_almost_equal(jd, jd_desired, decimal=12)


def test_thetaGMST_from_skyfield():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = pp.julian_date(dt)
    theta, thetadt = pp.theta_GMST1982(jd)
    assert_almost_equal(theta, 5.459562584754709, decimal=15)

# def test_sgp4():
#     """
#     Vallado, Eg.11-6 using SGP4 as propagator
#     """
#     tle1 = '1 16609U 86017A   93352.53502934  .00007889  0000  0  10529-3 0   342'
#     tle2 = '2 16609  51.6190  13.3340 0005770 102.5680 257.5950 15.59114070447869'
#     sat = pp.twoline2rv(tle1, tle2, pp.wgs72)
#     from pprint import pprint
#     pprint(sat)
#     r0T = np.array([6585.038266, 1568.184321, 9.116355])
#     v0T = np.array([-1.1157766, 4.6316816, 6.0149576])
#     jd0T = 2450540.400
#     assert_almost_equal(pp.julian_date(sat.epoch),jd0T)


def test_IJK2SEZ():
    """
    Curtis, Eg. 5.9, p.270
    """
    from math import sin, cos, radians, asin, acos, degrees
    # Satellite position
    r = np.array([-2032.4, 4591.2, -4544.8])  # km, geocentric equatorial pos.
    # Location position
    H = 0   # elevation, sea level
    phi = -40    # deg latitude
    theta = 110  # deg, local sidereal time
    R_obs = np.array([-1673.0, 4598.0, -4078.0])  # km, from Eq. 5.56

    # sez = pp.IJK2SEZ(r, phi, theta, H)
    rho = r - R_obs  # relative position vector
    phi_rad = radians(phi)
    theta_rad = radians(theta)

    # rotation matrix, from Eq. 5.62a
    Q = np.array([
       [-sin(theta_rad), cos(theta_rad), 0.],
       [-sin(phi_rad)*cos(theta_rad), -sin(phi_rad)*sin(theta_rad), cos(phi_rad)],
       [cos(phi_rad)*cos(theta_rad), cos(phi_rad)*sin(theta_rad), sin(phi_rad)]
    ])

    rho2 = np.dot(Q, rho)

    rho2_hat = (1/np.linalg.norm(rho))*rho2

    a = degrees(asin(rho2_hat[2]))
    A = degrees(acos(rho2_hat[1]/cos(radians(a))))

    A_true = 129.8  # deg, azimuth
    a_true = 41.41  # deg, angular elevation
    assert_almost_equal(A, A_true, decimal=1)
    assert_almost_equal(a, a_true, decimal=1)


# def test_IJK2SEZ_2():
#     """
#     Vallado, p. 164
#     """

#     # example from Vallado, p.913
#     phi = 42.38    # degrees latitude
#     lmda = -71.13  # degrees local sidereal time
#     h_ellp = 24

#     r_sat_ECEF = np.array([885.7296, -4389.3856, 5069.5763])  # km
#     r_site_ECEF = pp.site_ECEF(phi, lmda, h_ellp)

#     from math import radians, degrees, sin, cos
#     # Rotation matrix IJK -> SEZ
#     angle2 = radians(90-phi)
#     angle3 = radians(lmda)

#     R = np.matmul(pp.rot2(angle2), pp.rot3(angle3))

#     # Get position vector from site to satellite
#     rho_ECEF = r_sat_ECEF - r_site_ECEF

#     # Rotate position vector to SEZ coordinates
#     rho_SEZ = np.matmul(R, rho_ECEF)

#     rho_SEZ_exact = np.array([-773.8654, -581.4980, 328.8145])  # km

#     assert_allclose(rho_SEZ, rho_SEZ_exact)

if __name__ == "__main__":
    # test_julian_date_vectorized()
    # test_theta_GMST1982()
    # test_appendix_c_conversion_from_TEME_to_ITRF_UTC1()
    # test_jd_from_skyfield2()
    # test_jd_from_skyfield3()
    # test_site_declination_and_K()
    # test_site_ECEF2()
    # test_ECEF_to_SEZ()
#     test_satellite_visible()

    def azm(s, e):
        out = np.arctan2(s, e) * 180/np.pi + 90
        print(out)
        if s<0 and e<0:
            out = out % 360
        return out

    s, e = -1, 1
    print(azm(s, e))

