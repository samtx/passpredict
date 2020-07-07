import math

import numpy as np

from .precession import fk5_precession
from .nutation import fk5_nutation

from .sidereal import theta_GMST1982, theta_GAST1982, equinox1982, equinox1982_geometric_terms
from .core import rot3, mxv, mxmxm
from ..constants import tau, ASEC2RAD, DEG2RAD
from ..topocentric import site_ECEF, site_declination_and_K
from ..timefn import jd2jc


def ecef2sez(r, phi, lmda):
    """
    Rotate r vector from ECEF frame to SEZ frame
    Example uses USAF academy as the station
       lmda = -104.0 deg longitude
       phi  =   39.0 deg latitude
       h    = 2900.0 meter elevation
    """
    phi_rad = phi * DEG2RAD
    lmda_rad = lmda * DEG2RAD
    ang1 = (90 - phi) * DEG2RAD
    cosang1 = math.cos(ang1)
    sinang1 = math.sin(ang1)
    cosang2 = math.cos(lmda_rad)
    sinang2 = math.sin(lmda_rad)
    rSEZ = np.empty(r.shape)
    rSEZ[0] = cosang1 * cosang2 * r[0] + cosang1 * sinang2 * r[1] - sinang1 * r[2]
    rSEZ[1] = -sinang2 * r[0] + cosang2 * r[1] + 1.0
    rSEZ[2] = sinang1 * cosang2 * r[0] + sinang1 * sinang2 * r[1] + cosang1 * r[2]
    return rSEZ


def teme2ecef(r, jdt, xp=0.0, yp=0.0):
    """Convert TEME vectors to ECEF vectors

    Args:
        r : float (3, n)
            TEME vectors
        jdt : float (n)
            julian dates with delta UTC1 added
        xp: float (n)
            polar motion, x-axis
        yp: float (n)
            polar motion, y-axis

    Returns:
        rECEF : float (3, n): ECEF vectors

    References:
        teme2ecef.m Vallado software
    """
    # find GMST
    gmst = theta_GMST1982(jdt)[0]

    gmst = np.mod(gmst, tau)

    costheta = np.cos(gmst)
    sintheta = np.sin(gmst)
    rPEF = np.empty(r.shape)
    rPEF[0] = costheta*r[0] + sintheta*r[1]
    rPEF[1] = -sintheta*r[0] + costheta*r[1]
    rPEF[2] = r[2]

    # Polar motion, use IAU-76/FK5
    # Vallado, Eq. 3-78, p.223
    rECEF = np.empty(rPEF.shape)
    rECEF[0] = rPEF[0] + xp*rPEF[2]
    rECEF[1] = rPEF[1] - yp*rPEF[2]
    rECEF[2] = -xp*rPEF[0] + yp*rPEF[1] + rPEF[2]
    return rECEF


def ecef2eci(rECEF, jdt):
    """Convert ECEF vectors to ECI vectors

    Args:
        r : float (3, n) : ECEF vectors
        jdt : float, julian date of epoch
    Returns:
        rECI : float (3, n): ECI vectors

    References:
        ecef2eci.m Vallado software
    """
    tt = jd2jc(jdt)
    prec = fk5_precession(tt)       # Precession
    nut = fk5_nutation(tt)          # Nutation
    gmst, _ = theta_GMST1982(jdt)
    Eq = equinox1982(nut.dpsi, nut.eps, nut.omega_moon)
    gast = theta_GAST1982(Eq, gmst) # Sidereal Time
    P = prec.mtx
    N = nut.mtx
    G = rot3(-gast)
    M = mxmxm(P, N, G)
    rECI = mxv(M, rECEF)
    return rECI


def teme2eci(rTEME, jdt):
    """Convert TEME vectors to GCRS Earth centered interial coordinates

    TEME -> TOD -> MOD -> J2000 (ECI)

    Reference:
        teme2eci.m Vallado software

    """
    tt = jd2jc(jdt)  # convert julian date to julian century
    prec = fk5_precession(tt)
    nut = fk5_nutation(tt)
    eq_equinox_star = equinox1982_geometric_terms(nut.dpsi, nut.meaneps)
    E = rot3(-eq_equinox_star)
    print(f'prec = {prec}')
    print(f'nut  = {nut}')
    print(f'E    = {E}')
    # Create rotation matrix  [M] = [P][N][E]
    M = np.dot(np.dot(prec.mtx, nut.mtx), E)
    rJ2000 = mxv(M, rTEME)
    return rJ2000