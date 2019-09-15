import numpy as np
from numpy import dot, cross
from math import sqrt, sin, cos, cosh, acosh, tan, atan, acos, radians, degrees, pi
from datetime import datetime

# Constants
R_EARTH = np.float64(6378.137)  # [km] Mean equatorial radius
R2_EARTH = np.longdouble(40680631.59076899)  # [km2] mean equat. radius sq.
e_EARTH = np.longdouble(0.081816221456)  # Earth eccentricity
e2_EARTH = np.longdouble(0.006694385000)  # Earth eccentricity squared
MU = np.longdouble(398600.4418)  # [km3/(solar s)2] gravitational parameter
J2 = np.longdouble(0.0010826267)


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


def fk5():
    pass


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
                    [0., cos(a), 1.]])
    return mtx


def coe2rv(p, e, i, Omega, w, nu, u=0., lmda_true=0., w_hat_true=0.):
    """Compute ECI position and velocity vectors from classical orbital elements

    References:
        Vallado, Algorithm 10, p.118
    """
    # special cases
    small = 1e-10
    i_rad = radians(i)
    if e < small:
        # circular
        if (i_rad < small) or (np.abs(i_rad - pi) < small):
            # equatorial
            w, Omega, nu = 0., 0., lmda_true
        else:
            # inclined
            w, nu = 0., u
    else:
        # elliptical
        if (i_rad < small) or (np.abs(i_rad - pi) < small):
            # equatorial
            Omega, w = 0., w_hat_true

    nu_rad = radians(nu)
    cos_nu = cos(nu_rad)
    sin_nu = sin(nu_rad)
    sqrt_MU_p = sqrt(MU / p)

    rPQW = np.array(
        [p * cos_nu / (1 + e * cos_nu), p * sin_nu / (1 + e * cos_nu), 0.])
    vPQW = np.array([-sqrt_MU_p * sin_nu, sqrt_MU_p * (e + cos_nu), 0.])

    Omega_rad = radians(Omega)
    w_rad = radians(w)
    cos_Om = cos(Omega_rad)
    cos_w = cos(w_rad)
    cos_i = cos(i_rad)
    sin_Om = sin(Omega_rad)
    sin_w = sin(w_rad)
    sin_i = sin(i_rad)

    rot_mtx = np.array([[
        cos_Om * cos_w - sin_Om * sin_w * cos_i,
        -cos_Om * sin_w - sin_Om * cos_w * cos_i, sin_Om * sin_i
    ],
                        [
                            sin_Om * cos_w + cos_Om * sin_w * cos_i,
                            -sin_Om * sin_w + cos_Om * cos_w * cos_i,
                            -cos_Om * sin_i
                        ], [sin_w * sin_i, cos_w * sin_i, cos_i]])

    rIJK = dot(rot_mtx, rPQW)
    vIJK = dot(rot_mtx, vPQW)

    return rIJK, vIJK


def rv2coe(rIJK, vIJK, findall=False):
    """Compute classical orbital elements from ECI position and velocity vectors.

    References:
        Vallado, Algorithm 9, p.113
    """
    hvec = cross(rIJK, vIJK)
    h = sqrt(dot(hvec, hvec))
    Khat = np.array([0., 0., 1.])
    nvec = cross(Khat, hvec)
    rmag = sqrt(dot(rIJK, rIJK))
    v2 = dot(vIJK, vIJK)
    rmaginv = 1.0 / rmag
    evec = (1 / MU) * ((v2 - MU * rmaginv) * rIJK - dot(rIJK, vIJK) * vIJK)
    e = sqrt(dot(evec, evec))

    xi = v2 * 0.5 - MU * rmaginv

    if e != 1.0:
        a = -MU / (2 * xi)
        p = a * (1 - e**2)
    else:
        p = h**2 / MU
        p = np.inf

    cos_i = hvec[2] / h
    i = degrees(acos(cos_i))

    n = sqrt(dot(nvec, nvec))
    cos_Om = nvec[0] / n
    cos_w = dot(nvec, evec) / (n * e)
    cos_nu = dot(evec, rIJK) * rmaginv * (1 / e)

    if nvec[1] < 0.:
        Omega = 360 - degrees(acos(cos_Om))
    else:
        Omega = degrees(acos(cos_Om))

    if evec[2] < 0.:
        w = 360 - degrees(acos(cos_w))
    else:
        w = degrees(acos(cos_w))

    if dot(rIJK, vIJK) < 0.:
        nu = 360 - degrees(acos(cos_nu))
    else:
        nu = degrees(acos(cos_nu))

    # special cases
    small = 1e-10
    what_true, lmda_true, u = 0., 0., 0.
    # elliptical equatorial
    if ((i < small or abs(pi - i) < small) and (e >= small)) or findall:
        cos_what_true = evec[0] / e
        if evec[1] < 0:
            what_true = 360 - degrees(acos(cos_what_true))
        else:
            what_true = degrees(acos(cos_what_true))
    # circular inclined
    if ((i >= small or abs(pi - i) >= small) and (e < small)) or findall:
        cos_u = dot(nvec, rIJK) * (1 / n) * rmaginv
        if rIJK[2] < 0:
            u = 360 - degrees(acos(cos_u))
        else:
            u = degrees(acos(cos_u))
    # circular equatorial
    if ((i < small or abs(pi - i) < small) and (e < small)) or findall:
        cos_lmda_true = rIJK[0] * rmaginv
        if rIJK[1] < 0:
            lmda_true = 360 - degrees(acos(cos_lmda_true))
        else:
            lmda_true = degrees(acos(cos_lmda_true))

    # dictionary of classical orbital elements
    coe = {}
    coe['p'] = p
    coe['a'] = a
    coe['e'] = e
    coe['i'] = i
    coe['Omega'] = Omega
    coe['w'] = w
    coe['nu'] = nu
    coe['u'] = u
    coe['lmda_true'] = lmda_true
    coe['what_true'] = what_true

    return coe


def julian_date(yr, mo=None, dy=None, hr=None, mn=None, sec=None):
    """Compute Julian Date from datetime or equivalent elements

    References:
        Vallado, Algorithm 14, p.183
    """
    if isinstance(yr, datetime):
        dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
    jd1 = float(367 * yr)
    jd2 = float(int(7 * (yr + int((mo + 9) / 12)) / 4))
    jd3 = float(int(275 * mo / 9))
    jd4 = float(dy)
    jd5 = 1721013.5
    jd6 = float(((sec / 60 + mn) / 60 + hr) / 24)
    jd = jd1 - jd2 + jd3 + jd4 + jd5 + jd6
    print([jd1, jd2, jd3, jd4, jd5, jd6])
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


def pkepler_rv(r0, v0, dt, ndt=0., nddt=0.):
    """Propagate object from initial position and velocity vectors

    References:
        Vallado, Algorithm 65, p.691
    """

    coe = rv2coe(r0, v0)

    p0 = coe.get('p')
    a0 = coe.get('a')
    e0 = coe.get('e')
    i0 = coe.get('i')
    Omega0 = coe.get('Omega')
    w0 = coe.get('w')
    u = coe.get('u')
    lmda_true = coe.get('lmda_true')
    nu0 = coe.get('nu')

    if not a0:
        a0 = 1 / p0 * (1 - e0**2)

    if e0 != 0.0:
        E0 = nu2anomaly(nu0, e0)
    else:
        if u:
            E0 = u
        else:
            E0 = lmda_true

    M0 = E0 - e0 * sin(radians(E0))
    p0 = a0 * (1 - e0**2)
    n0 = sqrt(MU / a0**3)

    # Update for permutations
    # breakpoint()
    a = a0 - 2 * a0 / (3 * n0) * ndt * dt
    e = e0 - 2 * (1 - e0) / (3 * n0) * ndt * dt
    Omega = Omega0 - 3 * n0 * R2_EARTH * J2 / (2 * p0**2) * cos(
        radians(i0)) * dt
    w = w0 + 3 * n0 * R2_EARTH * J2 / (4 * p0**2) * (
        4 - 5 * sin(radians(i0))**2) * dt
    M = M0 + n0 * dt + ndt / 2. * dt**2 + nddt / 6. * dt**3
    p = a * (1 - e**2)

    E = kepEqtnE(M, e)
    if e != 0.0:
        nu = anomaly2nu(E, e)
        u, lmda_true = 0., 0.
    else:
        nu = nu0
        u, lmda_true = E, E

    r, v = coe2rv(p, e, i0, Omega, w, nu, u, lmda_true)

    return r, v


def pkepler_coe(a0,
                e0,
                i0,
                Omega0,
                w0,
                nu0,
                u,
                lmda_true,
                what_true,
                dt,
                ndt=0.,
                nddt=0.):
    """Propagate object from orbital elements using Kepler's equation

    References:
        Vallado, Algorithm 65, p.691
    """

    if e0 != 0.0:
        E0 = nu2anomaly(nu0, e0)
    else:
        if u:
            E0 = u
        else:
            E0 = lmda_true

    M0 = E0 - e0 * sin(radians(E0))
    p0 = a0 * (1 - e0**2)
    n0 = np.sqrt(MU / a0**3)

    # Update for permutations
    a = a0 - 2 * a0 / (3 * n0) * ndt * dt
    e = e0 - 2 * (1 - e0) / (3 * n0) * ndt * dt
    Omega = Omega0 - 3 * n0 * R2_EARTH * J2 / (2 * p0**2) * cos(
        radians(i0)) * dt
    w = w0 + 3 * n0 * R2_EARTH * J2 / (4 * p0**2) * (
        4 - 5 * sin(radians(i0))**2) * dt
    M = M0 + n0 * dt + ndt / 2. * dt**2 + nddt / 6. * dt**3
    p = a * (1 - e**2)

    E = kepEqtnE(M, e)
    if e != 0.0:
        nu = anomaly2nu(E, e)
        u, lmda_true = 0., 0.
    else:
        nu = nu0
        u, lmda_true = E, E

    r, v = coe2rv(p, e, i0, Omega, w, nu, u, lmda_true)
    return r, v


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


_second = 1.0 / (24.0 * 60.0 * 60.0)


def theta_GMST1982(jd_ut1):
    """Return the angle of Greenwich Mean Standard Time 1982 given the JD.

    This angle defines the difference between the idiosyncratic True
    Equator Mean Equinox (TEME) frame of reference used by SGP4 and the
    more standard Pseudo Earth Fixed (PEF) frame of reference.

    From AIAA 2006-6753 Appendix C.

    """
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
    angular_velocity = array([zero, zero, -theta_dot])
    R = rot_z(-theta)

    if len(rTEME.shape) == 1:
        rPEF = (R).dot(rTEME)
        vPEF = (R).dot(vTEME) + cross(angular_velocity, rPEF)
    else:
        rPEF = einsum('ij...,j...->i...', R, rTEME)
        vPEF = einsum('ij...,j...->i...', R, vTEME) + cross(
            angular_velocity, rPEF, 0, 0).T

    if xp == 0.0 and yp == 0.0:
        rITRF = rPEF
        vITRF = vPEF
    else:
        W = (rot_x(yp)).dot(rot_y(xp))
        rITRF = (W).dot(rPEF)
        vITRF = (W).dot(vPEF)
    return rITRF, vITRF


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
    import datetime
    import numpy as np
    from sgp4.propagation import sgp4 as sgp4_BR
    from sgp4.io import twoline2rv
    from sgp4.earth_gravity import wgs72, wgs84
    import matplotlib.pyplot as plt
    import pickle
    from skyfield.api import Topos, load
    from datetime import datetime, timedelta
    from pytz import timezone

    # Generate elevation function from Alfano using Skyfield
    recompute = False
    # sat = 'spirale a'
    sat = 'iss'
    sec_per_query = 1
    lofi_step = 60  # seconds
    dt_days = 2

    try:
        data = np.load('iss_razel.npz', allow_pickle=True)
        tt = data['tt']
        azs = data['azs']
        alts = data['alts']
        rr = data['rr']
        rrsat = data['rrsat']
    except:
        tend = datetime.datetime.now(
            datetime.timezone.utc) + datetime.timedelta(days=dt_days)
        (tt, azs, alts, rr, rrsat) = razel_from_t(tend,
                                                  sec_per_query=sec_per_query,
                                                  sat=sat)
        if isinstance(tt[0], datetime.datetime):
            # Convert datetime array to minutes since epoch (beginning)
            tmp = tt.copy()
            t1 = tmp[0]
            tt = np.empty(tmp.size)
            for i in range(tmp.size):
                dt = tmp[i] - t1
                tt[i] = dt.total_seconds() / 60
        np.savez('iss_razel.npz',
                 tt=tt,
                 azs=azs,
                 alts=alts,
                 rr=rr,
                 rrsat=rrsat)

    try:
        data = np.load('elev_fn.npz', allow_pickle=True)
        tt = data['tt']
        alts = data['alts']
        rr = data['rr']
        el_lower_lim = data['el_lower_lim']
        elev_fn = data['elev_fn']
        tt_lofi = data['tt_lofi']
        alts_lofi = data['alts_lofi']
    except:
        # Get elevation limit function from Alfano (1992)
        SEC_PER_HOUR = 3600
        start_hour = 0
        end_hour = 48
        tt = tt[SEC_PER_HOUR * start_hour:SEC_PER_HOUR * end_hour]
        alts = alts[SEC_PER_HOUR * start_hour:SEC_PER_HOUR * end_hour]
        tt_lofi = tt[::lofi_step]
        alts_lofi = alts[::lofi_step]
        el_lower_lim = 10   # deg lower limit
        elev_fn = alts_lofi - el_lower_lim
        np.savez('elev_fn.npz',
                 elev_fn=elev_fn,
                 tt=tt,
                 alts=alts,
                 rr=rr,
                 el_lower_lim=el_lower_lim,
                 tt_lofi=tt_lofi,
                 alts_lofi=alts_lofi)

    # Do the parabolic blending to find rise and set times
    C = lambda a, tau: a[3] * tau**3 + a[2] * tau**2 + a[1] * tau + a[0]
    T = np.linspace(0, 1)
    t = np.empty(0)
    el = np.empty(0)
    troot = np.empty(0)
    elroot = np.empty(0)
    from timeit import default_timer as tic
    t_0 = tic()
    for i in range(1, tt_lofi.size - 2):
        if (np.product(elev_fn[i-1:i+1]) > 0) and (np.product(elev_fn[i:i+2]) > 0) \
            and (np.product(elev_fn[i+1:i+3]) > 0):
            continue
        a_el = compute_alphas(elev_fn[i - 1], elev_fn[i], elev_fn[i + 1],
                              elev_fn[i + 2])
        el = np.append(el, C(a_el, T))
        a_t = compute_alphas(tt_lofi[i - 1], tt_lofi[i], tt_lofi[i + 1],
                             tt_lofi[i + 2])
        t = np.append(t, C(a_t, T))
        Troots = para_roots(a_el)
        # plot_para_blending()
        if not np.any(Troots):
            continue
        elroot = np.append(elroot, C(a_el, Troots))
        troot = np.append(troot, C(a_t, Troots))
    idx = (troot > tt[0]) & (troot < tt[-1])
    troot = troot[idx]
    elroot = elroot[idx]
    t_1 = tic()
    print(f'Time = {t_1-t_0:.6f} sec')

    # trise, tset = [], []
    # assert tt_lofi.size == alts_lofi.size
    # n = tt_lofi.size
    # for i in range(3,n-4):
    #     t, p = tt_lofi[i-1:i+3], alts_lofi[i-1:i+3]
    #     rise, set = parabolic_blending(t, p)
    #     if rise:
    #         trise.extend(rise)
    #     if set:
    #         tset.extend(set)
    # trise = np.asarray(trise)
    # tset = np.asarray(tset)

    ymax = alts.max()
    ymin = alts.min()

    fig, axs = plt.subplots(figsize=(30, 5), sharex=True)
    markersize = 1.25
    try:
        ax = axs[0]
    except:
        ax = axs
    ax.plot(tt, alts, ':', markersize=markersize, label='Exact Elev.')
    ax.plot(tt_lofi,
            elev_fn + el_lower_lim,
            'o-',
            markersize=markersize * 2,
            label='Elev. Fn.')
    ax.plot(tt_lofi,
            np.ones(tt_lofi.size) * el_lower_lim,
            'k--',
            label='Elev. Limit')
    ax.plot(t, el + el_lower_lim, 's-', label='Para Blend')
    ax.plot(troot,
            elroot + el_lower_lim,
            'o',
            markersize=markersize * 3,
            label='Para roots')
    ax.set_xlabel('Minutes since epoch')
    ax.set_ylabel('Elevation Angle')
    ax.grid()
    ax.legend()
    try:
        ax = axs[1]
        ax.plot(tt, rr, 'o-', markersize=markersize, label=r'$\rho$')
        ax.plot(tt, rrsat, 'o-', markersize=markersize, label=r'$R$')
        ax.legend()
        ax.grid()
    except:
        pass
    plt.savefig('iss_razel.png', dpi=400)
    with open('iss_razel.fig', 'wb') as figfile:
        pickle.dump(fig, figfile)
    plt.show()
