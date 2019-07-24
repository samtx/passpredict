import numpy as np
from numpy import dot, cross
from math import sqrt, sin, cos, asin, acos, radians, degrees, pi
from pprint import pprint

# Constants
R_EARTH = np.float64(6378.137)               # [km] Mean equatorial radius
R2_EARTH = np.longdouble(40680631.59076899)  # [km2] mean equat. radius sq.
e_EARTH = np.longdouble(0.081816221456)      # Earth eccentricity
e2_EARTH = np.longdouble(0.006694385000)     # Earth eccentricity squared
MU = np.longdouble(398600.4418)  # [km3/(solar s)2] gravitational parameter
J2 = np.longdouble(0.0010826267)


def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    phi_gd_rad = radians(np.float64(phi_gd))  # convert [deg] to [rad]
    h_ellp_km = np.float64(h_ellp)*(1/1000)   # convert [m] to [km]
    C = R_EARTH / np.sqrt(1 - e2_EARTH*sin(phi_gd_rad)**2)
    S = C*(1-e2_EARTH)
    r_delta = (C + h_ellp_km)*cos(phi_gd_rad)
    r_K = (S + h_ellp_km)*sin(phi_gd_rad)
    return (r_delta, r_K)


def site_ECEF(phi_gd, lmda, h_ellp):
    """
    Algorithm 51 from Vallado, p.430

    """

    r_delta, r_K = site_declination_and_K(phi_gd, h_ellp)

    lmda_rad = radians(lmda)

    r_site_ECEF = np.array([
        r_delta*cos(lmda_rad),
        r_delta*sin(lmda_rad),
        r_K
    ])

    return r_site_ECEF


def fk5():
    pass


def nu2anomaly(nu, e):
    """
    Vallado, p.77
    """
    nu_rad = radians(nu)
    if e < 1.0:
        # sinE = sin(nu_rad)*np.sqrt(1-e**2)/(1+e*cos(nu_rad))
        # E = asin(sinE)
        cosE = (e + cos(nu_rad))/(1 + e*cos(nu_rad))
        E = acos(cosE)
    return degrees(E)


def anomaly2nu(e, E):
    """
    Vallado, alg 6, p.77
    """
    E_rad = radians(E)
    if e < 1.0:
        cosNu = (cos(E_rad) - e)/(1 - e*cos(E_rad))
        nu = acos(cosNu)
    return degrees(nu)


def rot1(a):
    """
    Euler angle rotation matrix, first angle
    Vallado, Eq.3-15
    """
    mtx = np.array([
        [1.,      0.,     0.],
        [0.,  cos(a), sin(a)],
        [0., -sin(a), cos(a)]
    ])
    return mtx


def rot2(a):
    """
    Euler angle rotation matrix, second angle
    Vallado, Eq.3-15
    """
    mtx = np.array([
        [cos(a),     0., -sin(a)],
        [0.,     cos(a),      0.],
        [sin(a),     0.,  cos(a)]
    ])
    return mtx


def rot3(a):
    """
    Euler angle rotation matrix, third angle
    Vallado, Eq.3-15
    """
    mtx = np.array([
        [cos(a), sin(a),   0.],
        [-sin(a), cos(a),  0.],
        [0.,      cos(a),  1.]
    ])
    return mtx


def coe2rv(p, e, i, Omega, w, nu, u=0., lmda_true=0., w_hat_true=0.):
    """
    Vallado, Alg 10, p.118
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
    sqrt_MU_p = sqrt(MU/p)

    rPQW = np.array([
        p*cos_nu/(1 + e*cos_nu),
        p*sin_nu/(1 + e*cos_nu),
        0.
    ])
    vPQW = np.array([
        -sqrt_MU_p*sin_nu,
        sqrt_MU_p*(e + cos_nu),
        0.
    ])

    Omega_rad = radians(Omega)
    w_rad = radians(w)
    cos_Om = cos(Omega_rad)
    cos_w = cos(w_rad)
    cos_i = cos(i_rad)
    sin_Om = sin(Omega_rad)
    sin_w = sin(w_rad)
    sin_i = sin(i_rad)

    rot_mtx = np.array([
        [cos_Om*cos_w-sin_Om*sin_w*cos_i, -cos_Om*sin_w-sin_Om*cos_w*cos_i,  sin_Om*sin_i],
        [sin_Om*cos_w+cos_Om*sin_w*cos_i, -sin_Om*sin_w+cos_Om*cos_w*cos_i, -cos_Om*sin_i],
        [                    sin_w*sin_i,                      cos_w*sin_i,         cos_i]
    ])

    rIJK = dot(rot_mtx, rPQW)
    vIJK = dot(rot_mtx, vPQW)

    return rIJK, vIJK


def rv2coe(rIJK, vIJK, findall=False):
    """
    Vallado, Alg 9, p.113
    """
    hvec = cross(rIJK, vIJK)
    h = sqrt(dot(hvec, hvec))
    Khat = np.array([0., 0., 1.])
    nvec = cross(Khat, hvec)
    rmag = sqrt(dot(rIJK, rIJK))
    v2 = dot(vIJK, vIJK)
    rmaginv = 1.0/rmag
    evec = (1/MU)*((v2 - MU*rmaginv)*rIJK - dot(rIJK, vIJK)*vIJK)
    e = sqrt(dot(evec, evec))

    xi = v2*0.5 - MU*rmaginv

    if e != 1.0:
        a = -MU/(2*xi)
        p = a*(1 - e**2)
    else:
        p = h**2/MU
        p = np.inf

    cos_i = hvec[2]/h
    i = degrees(acos(cos_i))

    n = sqrt(dot(nvec, nvec))
    cos_Om = nvec[0]/n
    cos_w = dot(nvec, evec)/(n*e)
    cos_nu = dot(evec, rIJK)*rmaginv*(1/e)

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
        cos_what_true = evec[0]/e
        if evec[1] < 0:
            what_true = 360 - degrees(acos(cos_what_true))
        else:
            what_true = degrees(acos(cos_what_true))
    # circular inclined
    if ((i >= small or abs(pi - i) >= small) and (e < small)) or findall:
        cos_u = dot(nvec, rIJK)*(1/n)*rmaginv
        if rIJK[2] < 0:
            u = 360 - degrees(acos(cos_u))
        else:
            u = degrees(acos(cos_u))
    # circular equatorial
    if ((i < small or abs(pi - i) < small) and (e < small)) or findall:
        cos_lmda_true = rIJK[0]*rmaginv
        if rIJK[1] < 0:
            lmda_true = 360 - degrees(acos(cos_lmda_true))
        else:
            lmda_true = degrees(acos(cos_lmda_true))

    return (p, a, e, i, Omega, w, nu, u, lmda_true, what_true)


def kepEqtnE(M, e):
    """
    Vallado, Alg 2, p.65
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
        E1 = E0 + (M_rad - E0 + e*sin(E0))/(1 - e*cos(E0)) 
        if np.abs(E1 - E0) < atol:
            break
        E0 = E1
        k += 1
    return degrees(E1)

def pkepler_coe(a0, e0, i0, Omega0, w0, nu0, u, lmda_true, what_true, dt, ndt=0., nddt=0.):
    """
    Vallado, Alg 65, p.691

    Propagate object from orbital elements
    """

    if e0 != 0.0:
        E0 = nu2anomaly(nu0, e0)

    M0 = E0 - e0*sin(radians(E0))
    p0 = a0*(1 - e0**2)
    n0 = np.sqrt(MU/a0**3)

    # Update for permutations
    a = a0 - 2*a0/(3*n0)*ndt*dt
    e = e0 - 2*(1 - e0)/(3*n0)*ndt*dt
    Omega = Omega0 - 3*n0*R2_EARTH*J2/(2*p0**2)*cos(radians(i0))*dt
    w = w0 + 3*n0*R2_EARTH*J2/(4*p0**2)*(4 - 5*sin(radians(i0))**2)*dt
    M = M0 + n0*dt + ndt/2.*dt**2 + nddt/6.*dt**3
    p = a*(1-e**2)

    E = kepEqtnE(M, e)
    if e != 0.0:
        nu = anomaly2nu(E, e)
        u, lmda_true = 0., 0.
    else:
        nu = nu0
        u, lmda_true = E, E

    r, v = coe2rv(p, e, i0, Omega, w, nu, u, lmda_true)
    return r, v

def sgp4():
    pass
