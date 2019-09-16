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

def test_site_ECEF():
    """
    Vallado, Eg 7-1, p.431
    """
    phi_gd = 39.007     # [deg]
    lmda = -104.883     # [deg]
    alt = 2187.         # [m]
    r_ECEF = pp.site_ECEF(phi_gd, lmda, alt)
    r_ECEFtrue = np.array([-1275.1219, -4797.9890, 3994.2975])
    assert_allclose(r_ECEF, r_ECEFtrue)


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


def test_rv2coe():
    """
    Vallado, Eg.2-5, p.114
    """
    rIJK = np.array([6524.834, 6862.875, 6448.296]) # [km]
    vIJK = np.array([4.901327, 5.533756, -1.976341]) # [km/s]
    p_true, a_true, e_true = 11067.790, 36127.343, 0.832853
    i_true, Omega_true, w_true = 87.870, 227.898, 53.38
    nu_true, u_true, lmda_true_true, what_true_true = 92.335, 145.60549, 55.282587, 247.806
    coe = pp.rv2coe(rIJK, vIJK, findall=True)
    assert_almost_equal(coe['p'], p_true, decimal=2)
    assert_almost_equal(coe['a'], a_true, decimal=2)
    assert_almost_equal(coe['e'], e_true, decimal=6)
    assert_almost_equal(coe['i'], i_true, decimal=3)
    assert_almost_equal(coe['Omega'], Omega_true, decimal=3)
    assert_almost_equal(coe['w'], w_true, decimal=2)
    assert_almost_equal(coe['nu'], nu_true, decimal=3)
    assert_almost_equal(coe['u'], u_true, decimal=0)
    assert_almost_equal(coe['lmda_true'], lmda_true_true, decimal=3)
    assert_almost_equal(coe['what_true'], what_true_true, decimal=3)


def test_rv2coe_to_coe2rv():
    rIJK = np.array([6524.834, 6862.875, 6448.296])   # [km]
    vIJK = np.array([4.901327, 5.533756, -1.976341])  # [km/s]
    coe = pp.rv2coe(rIJK, vIJK, findall=True)
    p = coe.get('p')
    e = coe.get('e')
    i = coe.get('i')
    Omega = coe.get('Omega')
    w = coe.get('w')
    u = coe.get('u')
    lmda_true = coe.get('lmda_true')
    nu = coe.get('nu')
    what_true = coe.get('what_true')
    r_calc, v_calc = pp.coe2rv(p, e, i, Omega, w, nu, u, lmda_true, what_true)
    assert_allclose(r_calc, rIJK)
    assert_allclose(v_calc, vIJK)


def test_coe2rv_to_rv2coe():
    p = 11067.790   # [km]
    e = 0.83285
    i = 87.87       # [deg]
    Omega = 227.89  # [deg]
    w = 53.38       # [deg]
    nu = 92.335     # [deg]
    r, v = pp.coe2rv(p, e, i, Omega, w, nu)
    coeT = pp.rv2coe(r, v, findall=True)
    assert_almost_equal(p, coeT['p'], decimal=2)
    assert_almost_equal(e, coeT['e'], decimal=6)
    assert_almost_equal(i, coeT['i'], decimal=3)
    assert_almost_equal(Omega, coeT['Omega'], decimal=3)
    assert_almost_equal(w, coeT['w'], decimal=2)
    assert_almost_equal(nu, coeT['nu'], decimal=3)


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
#     E = pp.kepEqtnE(M, e)
#     nu = pp.anomaly2nu(e, E)
#     r_calc, v_calc = pp.coe2rv(p, e, i, Omega, w, nu)
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
#     # jdf = pp.julian_date(1997, 4, 2, 13, 8, 0)
#     jdf = 2450540.5472
#     dt = jdf - jd0
#     rf, vf = pp.pkepler_rv(r0, v0, dt, ndt)
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


# def test_riseset():
#     """
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
    assert_almost_equal(E_calc, E_true)


def test_julian_date():
    """
    Vallado, eg.3-4
    """
    yr, mo, dy = 1996, 10., 26.
    hr, mn, sec = 14., 20., 0.
    jd = pp.julian_date(yr, mo, dy, hr, mn, sec)
    jdT = 2450383.09722222
    assert_almost_equal(jd, jdT)


def test_julian_date_datetime():
    """
    Vallado, eg.3-4
    """
    yr, mo, dy = 1996, 10, 26
    hr, mn, sec = 14, 20, 0
    dt = datetime(yr, mo, dy, hr, mn, sec)
    jd = pp.julian_date(dt)
    jdT = 2450383.09722222
    assert_almost_equal(jd, jdT)


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
    theta *= 180/np.pi  # convert from radians to degrees
    assert_almost_equal(theta, 152.578787810)  # degrees
    # assert_almost_equal(thetadt, 152.578787810)  # degrees


def test_theta_GMST1982_2():
    """Compute the Greenwich Mean Sidereal Time

    References
        Vallado, Eg. 3-15, p.230
    """
    # April 6, 2004, 07:51:28.386 009UTC, not UTC1
    dt = datetime(2004, 4, 6, 7, 51, 28, 386009)
    jd = pp.julian_date(dt)
    theta, thetadt = pp.theta_GMST1982(jd)
    theta *= 180/np.pi
    assert_almost_equal(theta, 312.8098943, decimal=3)


def test_appendix_c_conversion_from_TEME_to_ITRF():
    """Test TEME to ITRF conversion

    References:
        Vallado et al., Revision 2
        Rhodes, Skyfield library, test_earth_satellites.py
    """
    rTEME = np.array([5094.18016210, 6127.64465950, 6380.34453270])
    vTEME = np.array([-4.746131487, 0.785818041, 5.531931288])
    vTEME = vTEME * 24.0 * 60.0 * 60.0  # km/s to km/day

    # Apr 6, 2004,  07:51:28.386 UTC
    dt = datetime(2004, 4, 6, 7, 51, 28, 386000)
    jd = pp.julian_date(dt)

    # Polar motion
    xp = -0.140682  # arcseconds
    yp = 0.333309   # arcseconds
    ARCSEC2RAD = np.pi/180 * (1/3600)
    xp *= ARCSEC2RAD
    yp *= ARCSEC2RAD
    xp = yp = 0.0
    rITRF, vITRF = pp.TEME_to_ITRF(jd, rTEME, vTEME, xp, yp)

    print(rITRF)
    assert_almost_equal(rITRF[0], -1033.47938300, decimal=1)
    assert_almost_equal(rITRF[1], 7901.29527540, decimal=1)
    assert_almost_equal(rITRF[2], 6380.35659580, decimal=1)


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
    dt = datetime(2004, 4, 6, 7, 51, 28, 386000)
    dt += timedelta(microseconds=-439961)
    jd = pp.julian_date(dt)

    # Polar motion
    xp = -0.140682  # arcseconds
    yp = 0.333309   # arcseconds
    ARCSEC2RAD = np.pi/180 * (1/3600.0)
    xp *= ARCSEC2RAD
    yp *= ARCSEC2RAD
    # xp = yp = 0.
    rITRF, vITRF = pp.TEME_to_ITRF(jd, rTEME, vTEME, xp, yp)

    print(rITRF)
    assert_almost_equal(rITRF[0], -1033.47938300, decimal=0)
    assert_almost_equal(rITRF[1], 7901.29527540, decimal=1)
    assert_almost_equal(rITRF[2], 6380.35659580, decimal=1)


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
    test_appendix_c_conversion_from_TEME_to_ITRF_UTC1()
