# precession.py
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .constants import ASEC2RAD

@dataclass
class PrecessionParams:
    zeta: float
    theta: float
    z: float
    mtx: np.ndarray


@lru_cache(maxsize=128)
def fk5_precession_angles(Td):
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


@lru_cache(maxsize=128)
def fk5_precession(tt: float):
    """
    Precession parameters and matrix for FK5 reduction with IAU 1976 precession theory

    Converts MOD -> J2000 (ECI)
    """
    zeta, theta, z = fk5_precession_angles(tt)
    coszeta = np.cos(zeta)
    sinzeta = np.sin(zeta)
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    cosz = np.cos(z)
    sinz = np.sin(z)
    P = np.array(
        (
            ( costheta*cosz*coszeta - sinz*sinzeta,  sinz*costheta*coszeta + sinzeta*cosz,  sintheta*coszeta),
            (-sinzeta*costheta*cosz - sinz*coszeta, -sinz*sinzeta*costheta + cosz*coszeta, -sintheta*sinzeta),
            (                       -sintheta*cosz,                        -sintheta*sinz,          costheta)
        ), dtype=np.float64
    )
    prec_params = PrecessionParams(
        zeta=zeta,
        z=z,
        theta=theta,
        mtx=P
    )
    return prec_params


def fk5_precess_rotation(r, zeta, theta, z):
    """Perform precession rotations on position vector r
    From IAU 76 Precession Theory. Rotation from MOD to GCRF.

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
