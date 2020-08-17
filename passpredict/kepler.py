from math import sqrt, sin, cos, cosh, acosh, tan, atan, acos, radians, degrees, pi
import datetime
from collections import namedtuple
from typing import NamedTuple

import numpy as np
from numpy import dot, cross

from .constants import R2_EARTH, MU, J2

# Model of classical orbital elements
class COE(NamedTuple):
    p: float
    a: float
    e: float
    i: float
    Omega: float
    w: float
    nu: float
    u: float
    lmda_true: float
    what_true: float


def coe2rv(p, e, i, Omega, w, nu, u=0.0, lmda_true=0.0, w_hat_true=0.0):
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
            w, Omega, nu = 0.0, 0.0, lmda_true
        else:
            # inclined
            w, nu = 0.0, u
    else:
        # elliptical
        if (i_rad < small) or (np.abs(i_rad - pi) < small):
            # equatorial
            Omega, w = 0.0, w_hat_true

    nu_rad = radians(nu)
    cos_nu = cos(nu_rad)
    sin_nu = sin(nu_rad)
    sqrt_MU_p = sqrt(MU / p)

    rPQW = np.array([p * cos_nu / (1 + e * cos_nu), p * sin_nu / (1 + e * cos_nu), 0.0])
    vPQW = np.array([-sqrt_MU_p * sin_nu, sqrt_MU_p * (e + cos_nu), 0.0])

    Omega_rad = radians(Omega)
    w_rad = radians(w)
    cos_Om = cos(Omega_rad)
    cos_w = cos(w_rad)
    cos_i = cos(i_rad)
    sin_Om = sin(Omega_rad)
    sin_w = sin(w_rad)
    sin_i = sin(i_rad)

    rot_mtx = np.array(
        [
            [
                cos_Om * cos_w - sin_Om * sin_w * cos_i,
                -cos_Om * sin_w - sin_Om * cos_w * cos_i,
                sin_Om * sin_i,
            ],
            [
                sin_Om * cos_w + cos_Om * sin_w * cos_i,
                -sin_Om * sin_w + cos_Om * cos_w * cos_i,
                -cos_Om * sin_i,
            ],
            [sin_w * sin_i, cos_w * sin_i, cos_i],
        ]
    )

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
    Khat = np.array([0.0, 0.0, 1.0])
    nvec = cross(Khat, hvec)
    rmag = sqrt(dot(rIJK, rIJK))
    v2 = dot(vIJK, vIJK)
    rmaginv = 1.0 / rmag
    evec = (1 / MU) * ((v2 - MU * rmaginv) * rIJK - dot(rIJK, vIJK) * vIJK)
    e = sqrt(dot(evec, evec))

    xi = v2 * 0.5 - MU * rmaginv

    if e != 1.0:
        a = -MU / (2 * xi)
        p = a * (1 - e ** 2)
    else:
        p = h ** 2 / MU
        p = np.inf

    cos_i = hvec[2] / h
    i = degrees(acos(cos_i))

    n = sqrt(dot(nvec, nvec))
    cos_Om = nvec[0] / n
    cos_w = dot(nvec, evec) / (n * e)
    cos_nu = dot(evec, rIJK) * rmaginv * (1 / e)

    if nvec[1] < 0.0:
        Omega = 360 - degrees(acos(cos_Om))
    else:
        Omega = degrees(acos(cos_Om))

    if evec[2] < 0.0:
        w = 360 - degrees(acos(cos_w))
    else:
        w = degrees(acos(cos_w))

    if dot(rIJK, vIJK) < 0.0:
        nu = 360 - degrees(acos(cos_nu))
    else:
        nu = degrees(acos(cos_nu))

    # special cases
    small = 1e-10
    what_true, lmda_true, u = 0.0, 0.0, 0.0
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

    # object of classical orbital elements
    return COE(p=p, a=a, e=e, i=i, Omega=Omega, w=w, nu=nu, u=u, lmda_true=lmda_true, what_true=what_true)


def pkepler_rv(r0, v0, dt, ndt=0.0, nddt=0.0):
    """Propagate object from initial position and velocity vectors

    References:
        Vallado, Algorithm 65, p.691
    """

    coe = rv2coe(r0, v0)

    p0 = coe.p
    a0 = coe.a
    e0 = coe.e
    i0 = coe.i
    Omega0 = coe.Omega
    w0 = coe.w
    u = coe.u
    lmda_true = coe.lmda_true
    nu0 = coe.nu

    if not a0:
        a0 = 1 / p0 * (1 - e0 ** 2)

    if e0 != 0.0:
        E0 = nu2anomaly(nu0, e0)
    else:
        if u:
            E0 = u
        else:
            E0 = lmda_true

    M0 = E0 - e0 * sin(radians(E0))
    p0 = a0 * (1 - e0 ** 2)
    n0 = sqrt(MU / a0 ** 3)

    # Update for permutations
    # breakpoint()
    a = a0 - 2 * a0 / (3 * n0) * ndt * dt
    e = e0 - 2 * (1 - e0) / (3 * n0) * ndt * dt
    Omega = Omega0 - 3 * n0 * R2_EARTH * J2 / (2 * p0 ** 2) * cos(radians(i0)) * dt
    w = (
        w0
        + 3 * n0 * R2_EARTH * J2 / (4 * p0 ** 2) * (4 - 5 * sin(radians(i0)) ** 2) * dt
    )
    M = M0 + n0 * dt + ndt / 2.0 * dt ** 2 + nddt / 6.0 * dt ** 3
    p = a * (1 - e ** 2)

    E = kepEqtnE(M, e)
    if e != 0.0:
        nu = anomaly2nu(E, e)
        u, lmda_true = 0.0, 0.0
    else:
        nu = nu0
        u, lmda_true = E, E

    r, v = coe2rv(p, e, i0, Omega, w, nu, u, lmda_true)

    return r, v


def pkepler_coe(
    a0, e0, i0, Omega0, w0, nu0, u, lmda_true, what_true, dt, ndt=0.0, nddt=0.0
):
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
    p0 = a0 * (1 - e0 ** 2)
    n0 = np.sqrt(MU / a0 ** 3)

    # Update for permutations
    a = a0 - 2 * a0 / (3 * n0) * ndt * dt
    e = e0 - 2 * (1 - e0) / (3 * n0) * ndt * dt
    Omega = Omega0 - 3 * n0 * R2_EARTH * J2 / (2 * p0 ** 2) * cos(radians(i0)) * dt
    w = (
        w0
        + 3 * n0 * R2_EARTH * J2 / (4 * p0 ** 2) * (4 - 5 * sin(radians(i0)) ** 2) * dt
    )
    M = M0 + n0 * dt + ndt / 2.0 * dt ** 2 + nddt / 6.0 * dt ** 3
    p = a * (1 - e ** 2)

    E = kepEqtnE(M, e)
    if e != 0.0:
        nu = anomaly2nu(E, e)
        u, lmda_true = 0.0, 0.0
    else:
        nu = nu0
        u, lmda_true = E, E

    r, v = coe2rv(p, e, i0, Omega, w, nu, u, lmda_true)
    return r, v


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
    E1 = 1.0
    k = 0
    while k < 1000:
        E1 = E0 + (M_rad - E0 + e * sin(E0)) / (1 - e * cos(E0))
        if np.abs(E1 - E0) < atol:
            break
        E0 = E1
        k += 1
    return degrees(E1)


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
