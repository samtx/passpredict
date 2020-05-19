import numpy as np
import math
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
from math import pi
from passpredictor.topocentric import site_ECEF, site_declination_and_K


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


def site2eci(lat, lon, h, jdt):
    """Compute site ECI vector for jdt array

    TO DO: write this function

    """
    rsiteECI = 0
    return rsiteECI


def fk5(r, xp=0.0, yp=0.0):
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


def eps1982(tt):
    """Compute the average eps 1982

    Args:
        tt : terrestial time

    References:
        Vallado, p. 225, Eq. 3-81
    """
    return 23.439291 * DEG2RAD + (-0.0130042 + (-1.64e-7 + 5.04e-7 * tt) * tt) * tt


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


def omega_moon(tt):
    """Nutation parameters for the moon

    Args:
        tt : terrestial time

    References:
        Vallado, p. 225, Eq. 3-82
    """
    r = 360 * 60
    omega_moon_deg = 125.04455501
    omega_moon_deg += (-5 * r - 134.1361851 + (0.0020756 + 2.139e-6 * tt) * tt) * tt
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
    eq += 0.00264 * ASEC2RAD * np.sin(omega_moon)
    eq += 0.000063 * np.sin(2 * omega_moon)
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
    i = np.array(
        [1, 9, 31, 2, 10, 32, 11, 33, 34, 12, 35, 13, 36, 38, 37], dtype=np.int
    )
    A = np.array(
        [
            -171996,
            -13187,
            -2274,
            2062,
            1426,
            712,
            -517,
            -386,
            -301,
            217,
            -158,
            129,
            123,
            63,
            63,
        ],
        dtype=np.float64,
    )
    B = np.array(
        [
            -174.2,
            -1.6,
            -0.2,
            0.2,
            -3.4,
            0.1,
            1.2,
            -0.4,
            0.0,
            -0.5,
            0.0,
            0.1,
            0.0,
            0.1,
            0.0,
        ]
    )
    C = np.array(
        [92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, -1, -70, -53, -33, -2]
    )
    D = np.array(
        [8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0]
    )
    an = np.array(
        [
            [0, 0, 0, 0, 1],
            [0, 0, 2, -2, 2],
            [0, 0, 2, 0, 2],
            [0, 0, 0, 0, 2],
            [0, 1, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 1, 2, -2, 2],
            [0, 0, 2, 0, 1],
        ]
    )


def fk5_nutation(tt):
    """
    Ref: Table D-6, p. 1043
    """
    api = an1 * M_moon + an2 * M_sun + an3 * Um_moon + an4 * D_sun + an5 * Omega_moon
    psi_tmp = np.sum(A + B * tt) * np.sin(api)
    eps_tmp = np.sum(C + D * tt) * np.cos(api)

    return dPsi1980, dEps1980


def fk5_precession(Td):
    """Calculate precession angle Theta

    Assume JD2000 epoch, T0 = 0

    Args
        Td: terrestial time

    References:
        Vallado, p. 227, Eq. 3-87
    """
    # Td = (jdt - J2000)/36525
    zeta = (2306.2181 + (0.30188 + 0.017998 * Td) * Td) * Td
    theta = (2004.3109 + (-0.42665 - 0.041833 * Td) * Td) * Td
    z = (2306.2181 + (1.09468 + 0.018203 * Td) * Td) * Td
    # Convert to radians
    zeta *= ASEC2RAD
    theta *= ASEC2RAD
    z *= ASEC2RAD

    return zeta, theta, z


def precess_rotation(r, zeta, theta, z):
    """Perform precession rotations on position vector r
    From IAU 80 Precession Theory. Rotation from MOD to GCRF.

    Args:
        r : (3, n), n is the number of observations
        zeta, theta, z (n): precession angles in radians

    References:
        Vallado, p. 228, Eq. 3-89
    """
    coszeta = np.cos(zeta)
    sinzeta = np.sin(zeta)
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    cosz = np.cos(z)
    sinz = np.sin(z)

    rGCRF = np.empty(r.shape, dtype=np.float)
    rGCRF[0] = (
        (costheta * cosz * coszeta - sinz * sinzeta) * r[0]
        + (sinz * costheta * coszeta + sinzeta * cosz) * r[1]
        + sintheta * coszeta * r[2]
    )
    rGCRF[1] = (
        (-sinzeta * costheta * cosz - sinz * coszeta) * r[0]
        + (-sinz * sinzeta * costheta + cosz * coszeta) * r[1]
        - sintheta * sinzeta * r[2]
    )
    rGCRF[2] = -sintheta * cosz * r[0] - sintheta * sinz * r[1] + costheta * r[2]
    return rGCRF


def teme2ecef(r, jdt):
    """Convert TEME vectors to ECEF vectors

    Args:
        r : float (3, n) : TEME vectors
        jdt : float (n) : julian dates
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
    rPEF[0] = costheta*r[0] - sintheta*r[1]
    rPEF[1] = sintheta*r[0] + costheta*r[1]
    rPEF[2] = r[2]
    rECEF = rPEF

    return rECEF


def ecef2eci(r, jdt):
    """Convert ECEF vectors to ECI vectors

    Args:
        r : float (3, n) : ECEF vectors
        jdt : float (n) : julian dates
    Returns:
        rECI : float (3, n): ECI vectors

    References:
        teme2ecef.m Vallado software
    """
    # find GMST
    gmst = theta_GMST1982(jdt)

    gmst = np.mod(gmst, tau)

    costheta = np.cos(gmstg)
    sintheta = np.sin(gmstg)
    rPEF = np.empty(r.shape)
    rPEF[0] = costheta*r[0] - sintheta*r[1]
    rPEF[1] = sintheta*r[0] + costheta*r[1]
    rPEF[2] = r[2]
    rECEF = rPEF

    return rECEF


def rot1(a):
    """Compute Euler angle rotation matrix, first angle

    References:
        Vallado, Eq. 3-15
    """
    mtx = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, math.cos(a), math.sin(a)],
            [0.0, -math.sin(a), math.cos(a)],
        ]
    )
    return mtx


def rot2(a):
    """Compute Euler angle rotation matrix, second angle

    References:
        Vallado, Eq. 3-15
    """
    mtx = np.array(
        [
            [math.cos(a), 0.0, -math.sin(a)],
            [0.0, math.cos(a), 0.0],
            [math.sin(a), 0.0, math.cos(a)],
        ]
    )
    return mtx


def rot3(a):
    """Compute Euler angle rotation matrix, third angle

    References:
        Vallado, Eq. 3-15
    """
    mtx = np.array(
        [
            [math.cos(a), math.sin(a), 0.0],
            [-math.sin(a), math.cos(a), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    return mtx


def site_sat_rotations(lat, lon, h, rsatECEF):
    rsiteECEF = site_ECEF(lat, lon, h)
    rho = rsatECEF - np.atleast_2d(rsiteECEF).T
    rSEZ = ecef2sez(rho, lat, lon)
    return rSEZ


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
        vPEF = np.dot(R, vTEME) + np.cross(angular_velocity, rPEF)
    else:
        rPEF = np.einsum("ij...,j...->i...", R, rTEME)
        vPEF = (
            np.einsum("ij...,j...->i...", R, vTEME)
            + np.cross(angular_velocity, rPEF, 0, 0).T
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