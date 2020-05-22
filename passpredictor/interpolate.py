# Polynomial blending and interpolation functions
from math import radians, degrees
import numpy as np
from skyfield.api import load, Topos

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

