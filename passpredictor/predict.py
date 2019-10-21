import numpy as np
from numpy import dot, cross
import math
from math import sqrt, sin, cos, cosh, acosh, tan, atan, acos, radians, degrees, pi
import datetime
from numpy.linalg import norm
from passpredictor.constants import (
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
from .rotations import rot1, rot2, rot3, theta_GMST1982


def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    phi_gd_rad = radians(np.float64(phi_gd))  # convert [deg] to [rad]
    h_ellp_km = np.float64(h_ellp) * (1 / 1000)  # convert [m] to [km]
    C = R_EARTH / np.sqrt(1 - e2_EARTH * sin(phi_gd_rad) ** 2)
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
    r_site_ECEF = np.array([r_delta * cos(lmda_rad), r_delta * sin(lmda_rad), r_K])
    return r_site_ECEF


def site_ECEF2(phi_gd, lmda, h_ellp):
    """Compute ECEF coordinates for tracking site on Earth

    Args:
        phi_gd: (float) geodetic latitutde of site in degrees
        lmda: (float) east longitude of site in degrees
        h_ellp: (float) height above the reference ellipse in meters

    References:
        Vallado, p. 428, Eq. 7-1
    """
    phi_gd_rad = phi_gd * DEG2RAD
    lmda_rad = lmda * DEG2RAD
    cosphi = math.cos(phi_gd_rad)
    sinphi = math.sin(phi_gd_rad)
    C = R_EARTH / math.sqrt(1 - e2_EARTH * (sinphi ** 2))
    S = C * (1 - e2_EARTH)
    h_ellp *= 0.001  # convert to km
    tmp = (C + h_ellp) * cosphi
    r_site_ECEF = np.array(
        [tmp * math.cos(lmda_rad), tmp * math.sin(lmda_rad), (S + h_ellp) * sinphi]
    )
    return r_site_ECEF


def rng_el(r):
    """Get range and elevation from SEZ vector"""
    rng = np.linalg.norm(r, axis=1)
    el = np.arcsin(r[2] / rng)
    el *= RAD2DEG
    return el, rng


def razel(r):
    """Get range, azimuth, and elevation from SEZ vector"""
    rng = np.linalg.norm(r, axis=0)
    el = np.arcsin(r[2] / rng) * RAD2DEG
    az = (np.arctan2(r[0], r[1]) + pi * 0.5) * RAD2DEG
    idx = np.all([r[0] < 0,r[1] < 0], axis=0)
    az[idx] %= 360
    return rng, az, el


def sun_pos(t):
    """Compute the Sun position vector

    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    """
    t_ut1 = (t - 2451545.0) / 36525
    t_tdb = t_ut1
    lmda_Msun = (280.4606184 + 36000.77005361 * t_tdb) % 360
    # M_sun = (357.5291092 + 35999.05034*t_tdb) % 360
    M_sun = (357.5277233 + 35999.05034 * t_tdb) % 360
    lmda_eclp = lmda_Msun + 1.914666471 * np.sin(M_sun * DEG2RAD)
    lmda_eclp += 0.019994643 * np.sin(2 * M_sun * DEG2RAD)
    r_sun_mag = 1.000140612 - 0.016708617 * np.cos(M_sun * DEG2RAD)
    r_sun_mag -= 0.000139589 * np.cos(2 * M_sun * DEG2RAD)
    eps = 23.439291 - 0.0130042 * t_tdb
    coslmda = np.cos(lmda_eclp * DEG2RAD)
    sinlmda = np.sin(lmda_eclp * DEG2RAD)
    coseps = np.cos(eps * DEG2RAD)
    sineps = np.sin(eps * DEG2RAD)
    r = np.empty((3, t.size), dtype=np.float)
    r[0] = r_sun_mag * coslmda
    r[1] = r_sun_mag * coseps * sinlmda
    r[2] = r_sun_mag * sineps * sinlmda
    r *= AU_KM
    return r


def vector_angle(r1, r2):
    """Compute the angle between two vectors

        r1 and r2 : (3, n), n is the number of observations
    """
    numerator = np.einsum("ij,ij->j", r1, r2)  # vectorized dot product
    denominator = norm(r1, axis=0) * norm(r2, axis=0)
    out = np.arccos(numerator / denominator) * RAD2DEG
    return out


def satellite_visible(rsatECI, rsiteECI, rho, jdt):
    """Determine visibility of satellite from Earth"""
    visible = np.zeros(jdt.size)
    # First, determine if satellite is above the horizon
    # find indecies where rho[2] > 0 --> idx
    vis_idx = np.nonzero(rho[2] > 0)[0]
    # Loop over indecies for times that satellite is over horizon
    for i in range(len(vis_idx)):
        idx = vis_idx[i]
        jdt_i = jdt[idx]
        rsatECI_i = rsatECI[:, idx]
        rsiteECI_i = rsiteECI[:, idx]
        rho_i = rho[:, idx]
        # Check if site is in daylight, compute dot product of sun position.
        rsun = sun_pos(jdt_i)
        if len(rsun.shape) > 1:
            rsun = rsun.flatten()
        # TODO: compute ECI vector for site position
        if np.dot(rsun, rsiteECI) > 0:
            # site is in daylight
            visible[idx] = 1
        else:
            # If nighttime, check if satellite is illuminated or in shadow
            if is_sat_illuminated(rsatECI_i, rsun):
                # satellite is illuminated, compute razel coordinates
                visible[idx] = 3
            else:
                # satellite is in shadow
                visible[idx] = 2
    return visible


def sun_sat_angle(rsat, rsun):
    """Compute the sun-satellite angle
    Args:
        rsat : satellite position vector in ECI coordinates
        rsun : sun position vector in ECI coordinates
    Output:
        angle in radians between the two vectors
    References:
        Vallado, p. 912, Alg. 74
    """
    sinzeta = norm(np.cross(rsun, rsat, axisa=0, axisb=0)) / (norm(rsun) * norm(rsat))
    return np.arcsin(sinzeta)


def sun_sat_orthogonal_distance(rsat, zeta):
    """
    Args:
        rsat : satellite position vector in ECI coordinates
        zeta : angle in radians between the satellite and sun vectors
    Output:
        distance from satellite to center of Earth orthogonal to sun vector
    """
    return norm(rsat) * np.cos(zeta - math.pi * 0.5)


def is_sat_illuminated(rsat, rsun):
    zeta = sun_sat_angle(rsat, rsun)
    dist = sun_sat_orthogonal_distance(rsat, zeta)
    return dist > R_EARTH


def azm(s, e):
    """Compute azimuth from topocentric horizon coordinates SEZ
    Args:
        s : south vector from SEZ coordinate
        e : east vector from SEZ coordinate
    Output:
        azimuth angle in radians with 0 at north.
    """
    out = np.arctan2(s, e) + pi * 0.5
    if s < 0 and e < 0:
        out = out % (2 * pi)
    return out


def elev(z, rhomag):
    """Compute elevation angle from topocentric horizon coordinates SEZ
    Args:
        z : Z vector from SEZ coordinate
        rhomag : magnitude of SEZ coordinate vector
    Output:
        elevation angle in radians with 0 at horizon, pi/2 straight up
    """
    return np.arcsin(z / rhomag)


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
        rPEF = np.einsum("ij...,j...->i...", R, rTEME)
        vPEF = (
            np.einsum("ij...,j...->i...", R, vTEME)
            + cross(angular_velocity, rPEF, 0, 0).T
        )
    if xp == 0.0 and yp == 0.0:
        rITRF = rPEF
        vITRF = vPEF
    else:
        W = (rot1(yp)).dot(rot2(xp))
        W = W.T
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
    elev = np.arctan(rhoZ / np.sqrt(rhoS ** 2 + rhoE ** 2))
    rng = np.sqrt(rhoS ** 2 + rhoE ** 2 + rhoZ ** 2)
    f_azm_lim = rhoE - rhoS * np.tan()


def razel_from_t(tend=None, sec_per_query=120, sat="iss"):
    """Use skyfield to calculate azimuth, elevation from satellite object"""
    sat_data = {
        "iss": {
            "url": "http://celestrak.com/NORAD/elements/stations.txt",
            "name": "ISS (ZARYA)",
        },
        "spirale a": {
            "url": "http://celestrak.com/NORAD/elements/military.txt",
            "name": "SPIRALE A",
        },
    }
    satellites = load.tle(sat_data[sat]["url"])
    iss = satellites[sat_data[sat]["name"]]
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
    myLocation = {"lat": "30.2672 N", "lon": "097.7431 W", "tz": "US/Central"}
    tz = timezone(myLocation["tz"])
    location = Topos(myLocation["lat"], myLocation["lon"])
    difference = iss - location
    tt = []
    alts = []
    azs = []
    rr = []
    rrsat = []
    for t in T:
        R = iss.at(t).position.km
        rrsat.append(np.sqrt(np.sum(R ** 2)))
        alt, az, r = difference.at(t).altaz()
        alt, az = (x.dms()[0] + x.dms()[1] / 60 for x in [alt, az])  # extract degrees
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


def get_rho_SEZ_from_t(dt_min, datetime_end, lat, lon, h=0.0):
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
        return ((-x) ** (1 / 3)) * -1
    else:
        return x ** (1 / 3)


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
        cosphi = -b / (2 * sqrt(-(a ** 3) / 27))
        sinphi = sqrt(1 - cosphi ** 2)
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
    return a[3] * T ** 3 + a[2] * T ** 2 + a[1] * T + a[0]


def plot_para_blending(alph_p, p):
    fig, ax = plt.subplots()
    xT = np.linspace(0, 1)
    yT = C(alph_p, xT)
    ax.plot(xT, yT, "-", label="C(T)")
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
    dT = 1.0 - 0.0
    # troots = Troots*dt/dT
    # Compute parabolic blending for times
    alph_t = compute_alphas(t[0], t[1], t[2], t[3])
    C = lambda a, T: a[3] * T ** 3 + a[2] * T ** 2 + a[1] * T + a[0]
    troots = C(alph_t, Troots)
    print(f"troots = {troots}")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    xT = np.linspace(0, 1)
    xt = np.linspace(t[0], t[3])
    yT = C(alph_p, xT)
    yTroots = C(alph_p, Troots)
    ax.plot(xT, yT, "-", label="C(T)")
    ax.plot(Troots, yTroots, "o", label="C(Troot)")
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


def alfano_approx(t, lat, lon, h=0.0):
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
    f_elev = (
        np.arccos(np.cos(elev_lim) / R_sat)
        - elev_lim
        - np.arccos(np.cos(elev) / R_sat)
        + elev_lim
    )
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
    F = a[4] * T ** 4 + a[3] * T ** 3 + a[2] * T ** 2 + a[1] * T + a[0]  #  0 <= T <= 4


if __name__ == "__main__":
    import pickle
    import matplotlib.pyplot as plt
