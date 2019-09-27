import numpy as np
from numpy import dot, cross
from math import sqrt, sin, cos, cosh, acosh, tan, atan, acos, radians, degrees, pi
import datetime

# Constants
R_EARTH = np.float64(6378.137)  # [km] Mean equatorial radius
R2_EARTH = np.longdouble(40680631.59076899)  # [km2] mean equat. radius sq.
e_EARTH = np.longdouble(0.081816221456)  # Earth eccentricity
e2_EARTH = np.longdouble(0.006694385000)  # Earth eccentricity squared
MU = np.longdouble(398600.4418)  # [km3/(solar s)2] gravitational parameter
J2 = np.longdouble(0.0010826267)
J2000 = 2451545.0

# Various constants required by Skyfield
AU_M = 149597870700             # per IAU 2012 Resolution B2
AU_KM = 149597870.700
ASEC360 = 1296000.0
DAY_S = 86400.0
# Angles.
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321
tau = 6.283185307179586476925287  # lower case, for symmetry with math.pi


def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    phi_gd_rad = radians(np.float64(phi_gd))  # convert [deg] to [rad]
    h_ellp_km = np.float64(h_ellp) * (1 / 1000)  # convert [m] to [km]
    C = R_EARTH / np.sqrt(1 - e2_EARTH * sin(phi_gd_rad)**2)
    S = C * (1 - e2_EARTH)
    r_delta = (C + h_ellp_km) * cos(phi_gd_rad)
    r_K = (S + h_ellp_km) * sin(phi_gd_rad)
    return (r_delta, r_K)


def site_ECEF(phi_gd, lmda, h_ellp):
    """Compute ECEF coordinates for tracking site on Earth

    References:
        Vallado, Algorithm 51, p.430
    """

    r_delta, r_K = site_declination_and_K(phi_gd, h_ellp)

    lmda_rad = radians(lmda)

    r_site_ECEF = np.array(
        [r_delta * cos(lmda_rad), r_delta * sin(lmda_rad), r_K])

    return r_site_ECEF


def fk5(r, xp=0., yp=0.):
    """IAU-76 / FK5 reductions for polar motion, nutation, precession

    Args
        r : position vector in ECI (ITRF) coordinates
        xp : polar motion along x axis in radians
        yp : polar motion along y axis in radians

    References:
        Vallado, Alg. 24, p.228
    """

    # polar motion using small angle approximations
    # Ref: Vallado, Eq 3-78
    rW = np.empty(r.shape)
    rW[0] = (1 + xp) * r[0]
    rW[1] = (1 - yp) * r[1]
    rW[2] = (1 - xp + yp) * r[2]

    # nutation, from IAU-1980 Theory of Nutation

    # We need to find the precession and nutation angles:
    # z, Theta, zeta, epsilon, deltaPsi1980, epsbar1980, thetaGAST1982

def utc2tt(UTC, deltaAT=37., deltaUT1=0.):
    """Compute terrestial time from UTC

    deltaAT is posted annualy in the Astronomical Almanac
        As of 2019 the offset is 37 seconds.
    deltaUT1 is posted daily by the US Navy

    Returns:
        TT : dynamic terrestial time in Julian Centuries

    References:
        Vallado, Section 3.5.4, Alg. 16
    """
    assert abs(deltaUT1) < 1
    UT1 = UTC + deltaUT1 / DAY_S
    # atomic time
    TAI = UTC + deltaAT / DAY_S
    # terrestial time in Julian Days
    JDtt = TAI + 32.184 / DAY_S
    # Convert JD to Julian Centuries
    TT = (JDtt - J2000) / 36525
    return TT


def eps1982(tt):
    """Compute the average eps 1982

    Args:
        tt : terrestial time

    References:
        Vallado, p. 225, Eq. 3-81
    """
    return 23.439291*DEG2RAD + (-0.0130042 + (-1.64e-7 + 5.04e-7*tt) * tt) * tt


def omega_moon(tt):
    """Nutation parameters for the moon

    Args:
        tt : terrestial time

    References:
        Vallado, p. 225, Eq. 3-82
    """
    r = 360 * 60
    omega_moon_deg = 125.04455501
    omega_moon_deg += (-5*r - 134.1361851 + (0.0020756 + 2.139e-6*tt)*tt)*tt
    return omega_moon_deg * DEG2RAD


def equinox1982(dPsi1980, eps1980, omega_moon):
    """Return the equation of equinoxes for sidereal time
    for IAU-76/FK5 reductions

    Args:
        t : time

    References:
        Vallado, p.224
    """
    eq = dPsi1980 * np.cos(eps1980)
    eq += 0.00264*ASEC2RAD*np.sin(omega_moon)
    eq += 0.000063*np.sin(2*omega_moon)
    return eq


def theta_GAST1982(Eq, gmst):
    """Return the Greenwich apparent sidereal time

    Sidereal time for IAU-76/FK5 reductions

    Args:
        Eq : equinox1982
        gmst : Greenwich mean sidreal time

    References:
        Vallado, p.224
    """
    return Eq + gmst


def nutation_coeff():
    i = np.array([      1,      9,    31,    2,   10,  32,   11,   33,   34,   12, 35, 13, 36, 38, 37], dtype=np.int)
    A = np.array([-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301,  217, -158, 129, 123, 63, 63], dtype=np.float64)
    B = np.array([ -174.2,   -1.6,  -0.2,  0.2, -3.4, 0.1,  1.2, -0.4,  0.0, -0.5, 0.0, 0.1, 0.0, 0.1, 0.0 ])
    C = np.array([  92025,   5736,   977, -895,   54,  -7,  224,  200,  129,  -95, -1, -70, -53, -33, -2])
    D = np.array([    8.9,   -3.1,  -0.5,  0.5, -0.1, 0.0, -0.6,  0.0, -0.1,  0.3, 0.0, 0.0, 0.0, 0.0, 0.0])
    an = np.array([
        [ 0,  0,  0,  0,  1],
        [ 0,  0,  2, -2,  2],
        [ 0,  0,  2,  0,  2],
        [ 0,  0,  0,  0,  2],
        [ 0,  1,  0,  0,  0],
        [ 1,  0,  0,  0,  0],
        [ 0,  1,  2, -2,  2],
        [ 0,  0,  2,  0,  1]
    ])

def fk5_nutation(tt):
    """
    Ref: Table D-6, p. 1043
    """
    api = an1*M_moon + an2*M_sun + an3*Um_moon + an4*D_sun + an5*Omega_moon
    psi_tmp = np.sum(A + B*tt)*np.sin(api)
    eps_tmp = np.sum(C + D*tt)*np.cos(api)

    return dPsi1980, dEps1980

def fk5_precession(jdt):
    """Calculate precession angle Theta

    Assume JD2000 epoch, T0 = 0

    References:
        Vallado, p. 227, Eq. 3-87
    """
    Td = (jdt - J2000)/36525
    zeta = (2306.2181 + (0.30188 + 0.017998*Td)*Td)*Td
    Theta = (2004.3109 + (-0.42665 - 0.041833*Td)*Td)*Td
    z = (2306.2181 + (1.09468 + 0.018203*Td)*Td)*Td
    return zeta, Theta, z


def nu2anomaly(nu, e):
    """Compute anomaly from nu and eccentricity

    References:
        Vallado, p.77
    """
    if e < 1.0:
        cosNu = cos(radians(nu))
        cosE = (e + cosNu) / (1 + e * cosNu)
        E = degrees(acos(cosE))
        out = E
    elif float(e) == 1.0:
        B = tan(nu / 2)
        out = B
    else:
        cosNu = cos(radians(nu))
        coshH = (e + cosNu) / (1 + e * cosNu)
        H = degrees(acosh(coshH))
        out = H
    return out


def anomaly2nu(e, E, B=0, p=0, r=1, H=0):
    """Compute nu from anomaly and eccentricity.

    References:
        Vallado, Algorithm 6, p.77
    """
    if e < 1.0:
        cosE = cos(radians(E))
        cosNu = (cosE - e) / (1 - e * cosE)
    elif float(e) == 1.0:
        cosNu = (p - r) / r
    else:
        coshH = cosh(radians(H))
        cosNu = (coshH - e) / (1 - e * coshH)
    nu = acos(cosNu)
    return degrees(nu)


def rot1(a):
    """Compute Euler angle rotation matrix, first angle

    References:
        Vallado, Eq. 3-15
    """
    mtx = np.array([[1., 0., 0.], [0., cos(a), sin(a)], [0., -sin(a), cos(a)]])
    return mtx


def rot2(a):
    """Compute Euler angle rotation matrix, second angle

    References:
        Vallado, Eq. 3-15
    """
    mtx = np.array([[cos(a), 0., -sin(a)], [0., cos(a), 0.],
                    [sin(a), 0., cos(a)]])
    return mtx


def rot3(a):
    """Compute Euler angle rotation matrix, third angle

    References:
        Vallado, Eq. 3-15
    """
    mtx = np.array([[cos(a), sin(a), 0.], [-sin(a), cos(a), 0.],
                    [0., 0., 1.]])
    return mtx


def julian_date(yr, mo=None, dy=None, hr=None, mn=None, sec=None):
    """Compute Julian Date from datetime or equivalent elements

    Notes
        T0 = 2451545.0

    References:
        Vallado, Algorithm 14, p.183
    """
    if isinstance(yr, datetime.datetime) or isinstance(yr, np.datetime64):
        if isinstance(yr, np.datetime64):
            dt = yr.astype(datetime.datetime)
        else:
            dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
        sec += dt.microsecond * (10**-6)
    jd1 = 367 * yr
    jd2 = 7 * (yr + (mo + 9) // 12) // 4
    jd3 = (275 * mo) // 9
    jd4 = dy
    jd5 = 1721013.5
    jd6 = ((sec / 60 + mn) / 60 + hr) / 24
    jd = jd1 - jd2 + jd3 + jd4 + jd5 + jd6
    # print([jd1, jd2, jd3, jd4, jd5, jd6])
    return jd


def kepEqtnE(M, e):
    """

    References:
        Vallado, Algorithm 2, p.65
    """
    M_rad = radians(M)
    if (-pi < M_rad < 0) or (M_rad > pi):
        E0 = M_rad - e
    else:
        E0 = M_rad + e
    atol = 1e-8
    E1 = 1.
    k = 0
    while k < 1000:
        E1 = E0 + (M_rad - E0 + e * sin(E0)) / (1 - e * cos(E0))
        if np.abs(E1 - E0) < atol:
            break
        E0 = E1
        k += 1
    return degrees(E1)


def sgp4(tle1, tle2, t):
    """
    t: minutes since epoch
    """
    satrec = twoline2rv(tle1, tle2, wgs72)
    n = len(t)
    r = np.zeros((n, 3))  # position vector
    v = np.zeros((n, 3))  # velocity vector
    for i, ti in enumerate(t):
        r[i], v[i] = sgp4_BR(satrec, ti)
    return r, v


##################
# From sgp4lib.py in skyfield
###################

def ITRF_position_velocity_error(t):
    """Return the ITRF position, velocity, and error at time `t`.

    The position is an x,y,z vector measured in au, the velocity is
    an x,y,z vector measured in au/day, and the error is a vector of
    possible error messages for the time or vector of times `t`.

    """
    rTEME, vTEME, error = sgp4(tle, t)
    rTEME /= AU_KM
    vTEME /= AU_KM
    vTEME *= DAY_S
    rITRF, vITRF = TEME_to_ITRF(t.ut1, rTEME, vTEME)
    return rITRF, vITRF, error


def _at(self, t):
    """Compute this satellite's GCRS position and velocity at time `t`."""
    rITRF, vITRF, error = self.ITRF_position_velocity_error(t)
    rGCRS, vGCRS = ITRF_to_GCRS2(t, rITRF, vITRF)
    return rGCRS, vGCRS, rGCRS, error


def theta_GMST1982(jd_ut1):
    """Return the angle of Greenwich Mean Standard Time 1982 given the JD.

    This angle defines the difference between the idiosyncratic True
    Equator Mean Equinox (TEME) frame of reference used by SGP4 and the
    more standard Pseudo Earth Fixed (PEF) frame of reference.

    Args:
        jd_ut1 (float): Julian date since Jan 1, 2000

    Returns:
        theta (float): GMST in radians
        thetadt (float): time derivative of GMST in rad/s

    References:
        Vallado, et al. "Revisiting Spacetrack Report #3", AIAA, 2006-6753, Appendix C.
        Rhodes, Skyfield library, github.com/skyfielders/python-skyfield --> sgp4lib.py
    """
    tau = 2 * pi
    _second = 1.0 / (24.0 * 60.0 * 60.0)
    T0 = 2451545.0  # JD for Jan 1, 2000
    t = (jd_ut1 - T0) / 36525.0
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t
    dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t
    theta = (jd_ut1 % 1.0 + g * _second % 1.0) * tau
    theta_dot = (1.0 + dg * _second / 36525.0) * tau
    return theta, theta_dot


def TEME_to_ITRF(jd_ut1, rTEME, vTEME, xp=0.0, yp=0.0):
    """Convert TEME position and velocity into standard ITRS coordinates.

    This converts a position and velocity vector in the idiosyncratic
    True Equator Mean Equinox (TEME) frame of reference used by the SGP4
    theory into vectors into the more standard ITRS frame of reference.
    The velocity should be provided in units per day, not per second.

    From AIAA 2006-6753 Appendix C.
    """
    theta, theta_dot = theta_GMST1982(jd_ut1)
    zero = theta_dot * 0.0
    angular_velocity = np.array([zero, zero, -theta_dot])
    R = rot3(-theta).T

    if len(rTEME.shape) == 1:
        rPEF = np.dot(R, rTEME)
        vPEF = np.dot(R, vTEME) + cross(angular_velocity, rPEF)
    else:
        rPEF = np.einsum('ij...,j...->i...', R, rTEME)
        vPEF = np.einsum('ij...,j...->i...', R, vTEME) + cross(
            angular_velocity, rPEF, 0, 0).T

    if xp == 0.0 and yp == 0.0:
        rITRF = rPEF
        vITRF = vPEF
    else:
        W = (rot1(yp)).dot(rot2(xp))
        W = W.T
        rITRF = (W).dot(rPEF)
        vITRF = (W).dot(vPEF)
    return rITRF, vITRF


def julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day int."""
    janfeb = month < 3
    return (day
            + 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            - 32075)


def julian_date2(yr, mo=1, dy=1, hr=0, mn=0, sec=0.0):
    """Given a proleptic Gregorian calendar date, return a Julian date float."""
    if isinstance(yr, datetime.datetime) or isinstance(yr, np.datetime64):
        if isinstance(yr, np.datetime64):
            dt = yr.astype(datetime.datetime)
        else:
            dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
        sec += dt.microsecond * (10**-6)

    return julian_day(yr, mo, dy) - 0.5 + (
        sec + mn * 60.0 + hr * 3600.0) / DAY_S


def alfano_approx2():
    """
    From Alfano (1993), quartic blending

    Use 3 min time steps

    Example uses USAF academy as the station
       lmda = -104.0 deg longitude
       phi  =   39.0 deg latitude
       h    = 2900.0 meter elevation
    """

    rhoIJK = rsatIJK - rsiteIJK
    phi = 10.0  # degree geodetic latitude
    lmda = 100  # degree east longitude

    angle2 = radians(90 - phi)
    angle3 = radians(lmda)
    R = np.dot(rot2(angle2), rot3(angle3))
    rhoSEZ = np.dot(R, rhoIJK)
    rhoS, rhoE, rhoZ = rhoSEZ[0], rhoSEZ[1], rhoSEZ[2]
    # compute azimuth, elevation, and range
    azm = np.arctan(rhoE / -rhoS)
    elev = np.arctan(rhoZ / np.sqrt(rhoS**2 + rhoE**2))
    rng = np.sqrt(rhoS**2 + rhoE**2 + rhoZ**2)

    f_azm_lim = rhoE - rhoS * np.tan()


def razel_from_t(tend=None, sec_per_query=120, sat='iss'):
    """Use skyfield to calculate azimuth, elevation from satellite object"""
    sat_data = {
        'iss': {
            'url': 'http://celestrak.com/NORAD/elements/stations.txt',
            'name': 'ISS (ZARYA)'
        },
        'spirale a': {
            'url': 'http://celestrak.com/NORAD/elements/military.txt',
            'name': 'SPIRALE A'
        }
    }
    satellites = load.tle(sat_data[sat]['url'])
    iss = satellites[sat_data[sat]['name']]

    tstart = iss.epoch.utc_datetime()
    if not tend:
        tend = tstart + timedelta(hours=72)
    yr = int(tstart.year)
    mo = int(tstart.month)
    dy = int(tstart.day)
    hr = int(tstart.hour)
    mn = int(tstart.minute)
    sc = int(tstart.second)
    sc_interval = int((tend - tstart).total_seconds())
    ts = load.timescale()
    T = ts.utc(yr, mo, dy, hr, mn, range(sc, sc_interval, sec_per_query))

    # Find overpasses
    myLocation = {'lat': '30.2672 N', 'lon': '097.7431 W', 'tz': 'US/Central'}

    tz = timezone(myLocation['tz'])
    location = Topos(myLocation['lat'], myLocation['lon'])
    difference = iss - location

    tt = []
    alts = []
    azs = []
    rr = []
    rrsat = []
    for t in T:
        R = iss.at(t).position.km
        rrsat.append(np.sqrt(np.sum(R**2)))
        alt, az, r = difference.at(t).altaz()
        alt, az = (x.dms()[0] + x.dms()[1] / 60
                   for x in [alt, az])  # extract degrees
        t = t.utc_datetime()
        tt.append(t)
        alts.append(alt)
        azs.append(az)
        rr.append(r.km)
    tt = np.asarray(tt)
    azs = np.asarray(azs)
    alts = np.asarray(alts)
    rr = np.asarray(rr)
    rrsat = np.asarray(rrsat)

    return tt, azs, alts, rr, rrsat


def get_rho_SEZ_from_t(dt_min, datetime_end, lat, lon, h=0.):
    pass


def compute_alphas(p1, p2, p3, p4):
    """For cubic polynomials"""
    a0 = p2
    a1 = -0.5 * p1 + 0.5 * p3
    a2 = p1 - 2.5 * p2 + 2 * p3 - 0.5 * p4
    a3 = -0.5 * p1 + 1.5 * p2 - 1.5 * p3 + 0.5 * p4
    return np.asarray([a0, a1, a2, a3])


def realcbrt(x):
    """Return real result of cube root"""
    if x < 0:
        return ((-x)**(1 / 3)) * -1
    else:
        return x**(1 / 3)


def para_roots(alph_p):
    from math import sqrt, atan2, cos, degrees, radians, pi
    # Find roots T1, T2, T3 of C(T) in interval 0 <= T < 1
    # Rearrange coefficients
    P = alph_p[2] / alph_p[3]
    Q = alph_p[1] / alph_p[3]
    R = alph_p[0] / alph_p[3]
    a = 1 / 3 * (3 * Q - P * P)
    b = 1 / 27 * (2 * P * P * P - 9 * P * Q + 27 * R)
    delta = a * a * a / 27 + b * b / 4
    tol = 1e-8
    Troots = None

    # plot_para_blending(alph_p)

    if delta > tol:
        # Use Cardan's solution, Eq. C-31
        sqrt_delta = sqrt(delta)
        halfb = b / 2
        x1 = realcbrt(-halfb + sqrt_delta) + realcbrt(-halfb - sqrt_delta)
        if (x1 > 0) and (x1 < 1):
            Troots = np.asarray([x1])
    elif np.abs(delta) <= tol:
        x1 = 2 * realcbrt(-b / 2)
        x2 = x3 = realcbrt(b / 2)
        Troots = np.asarray([x1, x2, x3])
        Troots = Troots[(Troots > 0) & (Troots < 1)]
        if Troots.size == 0:
            return None
    elif delta < -tol:
        E0 = 2 * sqrt(-a / 3)
        cosphi = -b / (2 * sqrt(-(a**3) / 27))
        sinphi = sqrt(1 - cosphi**2)
        phi = atan2(sinphi, cosphi)
        z1 = E0 * cos(phi / 3)
        z2 = E0 * cos(phi / 3 + 2 * pi / 3)
        z3 = E0 * cos(phi / 3 + 4 * pi / 3)
        Pdiv3 = P / 3
        x1 = z1 - Pdiv3
        x2 = z2 - Pdiv3
        x3 = z3 - Pdiv3
        Troots = np.array([x1, x2, x3])
        Troots = Troots[(Troots > 0) & (Troots < 1)]
        if Troots.size == 0:
            return None
    else:
        return None
    return Troots


def C(a, T):
    return a[3] * T**3 + a[2] * T**2 + a[1] * T + a[0]


def plot_para_blending(alph_p, p):
    fig, ax = plt.subplots()
    xT = np.linspace(0, 1)
    yT = C(alph_p, xT)
    ax.plot(xT, yT, '-', label='C(T)')
    ax.grid()
    ax.legend()
    plt.show()


def parabolic_blending(t, p):
    """Parabolic blending approximation from Alfano (1992)
    pts = (t1, p1), ... , (t4, p4)

    References:
        Vallado, Appendix C.5.1 for finding cubic roots
    """
    from math import sqrt, atan2, cos, degrees, radians, pi
    alph_p = compute_alphas(p[0], p[1], p[2], p[3])
    Troots = para_roots(alph_p)
    if not Troots:
        return None
    # our roots are x1, x2, x3
    # rescale time dimension from (0,1) -> (t0, t3)
    dt = t[3] - t[0]
    dT = 1. - 0.
    # troots = Troots*dt/dT
    # Compute parabolic blending for times
    alph_t = compute_alphas(t[0], t[1], t[2], t[3])
    C = lambda a, T: a[3] * T**3 + a[2] * T**2 + a[1] * T + a[0]
    troots = C(alph_t, Troots)
    print(f'troots = {troots}')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    xT = np.linspace(0, 1)
    xt = np.linspace(t[0], t[3])
    yT = C(alph_p, xT)
    yTroots = C(alph_p, Troots)
    ax.plot(xT, yT, '-', label='C(T)')
    ax.plot(Troots, yTroots, 'o', label='C(Troot)')
    ax.grid()
    ax.legend()
    plt.show()
    # Check if it is a rise or set time
    trise = []
    tset = []
    for i in range(Troots.size):
        T = Troots[i]
        res = C(alph_p, T + 0.05)
        if res > 0:
            # troots[i] is a rise time
            trise.append(troots[i])
        else:
            # troots[i] is a rise time
            tset.append(troots[i])
    return trise, tset


def alfano_approx(t, lat, lon, h=0.):
    """
    From Alfano (1993), described in Vallado p. 914
    """
    elev_lim = radians(10)  # elevation lower limit

    rho = get_rho_SEZ_from_t(dt_min, datetime_end, lat=lat, lon=lon, h=h)
    elev = get_elev_from_rho(rho)
    R_sat = get_Rsat_norm(dt_min, datetime_end)
    # def f_range(i):
    #     return sqrt(rho[i,0]**2+rho[i,1]**2+rho[i,2]**2) - rho_lim

    # def f_beta_lim(i):
    #     return rho[i,1]

    # Lower limit on elevation for viewing satellite, vectorized
    f_elev = np.arccos(np.cos(elev_lim) / R_sat) - elev_lim - np.arccos(
        np.cos(elev) / R_sat) + elev_lim

    # Find indecies where f_elev changes sign
    # from https://stackoverflow.com/questions/2652368/how-to-detect-a-sign-change-for-elements-in-a-numpy-array
    idx = np.where(np.sign(f_elev[:-1]) != np.sign(f_elev[1:]))[0] + 1

    # Note: could also use scipy.interpolate.UnivariateSpline(k=4 or 5) to create quartic and quintic polynomials and splines
    a = np.zeros(5)
    a[0] = p[0]
    a[1] = (-50 * p[0] + 96 * p[1] - 72 * p[2] + 32 * p[3] - 6 * p[4]) / 24
    a[2] = (35 * p[0] - 104 * p[1] + 114 * p[2] - 56 * p[3] + 11 * p[4]) / 24
    a[3] = (-10 * p[0] + 36 * p[1] - 48 * p[2] + 28 * p[3] - 6 * p[4]) / 24
    a[4] = (p[0] - 4 * p[1] + 6 * p[2] - 4 * p[3] + p[4]) / 24

    # scale t from t1 < t5 to 0 <= T <= 4
    F = a[4] * T**4 + a[3] * T**3 + a[2] * T**2 + a[1] * T + a[0]  #  0 <= T <= 4


if __name__ == "__main__":
    import pickle
    import matplotlib.pyplot as plt
