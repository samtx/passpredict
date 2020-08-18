import math

import numpy as np

from .constants import DEG2RAD, R_EARTH, e2_EARTH


def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    phi_gd_rad = phi_gd * DEG2RAD  # convert [deg] to [rad]
    h_ellp_km = h_ellp * 0.001  # convert [m] to [km]
    C = R_EARTH / math.sqrt(1 - e2_EARTH * math.sin(phi_gd_rad) ** 2)
    S = C * (1 - e2_EARTH)
    r_delta = (C + h_ellp_km) * math.cos(phi_gd_rad)
    r_K = (S + h_ellp_km) * math.sin(phi_gd_rad)
    return (r_delta, r_K)


def site_ECEF(phi_gd, lmda, h_ellp):
    """Compute ECEF coordinates for tracking site on Earth

    References:
        Vallado, Algorithm 51, p.430
    """
    r_delta, r_K = site_declination_and_K(phi_gd, h_ellp)
    lmda_rad = math.radians(lmda)
    r_site_ECEF = np.array([r_delta * math.cos(lmda_rad), r_delta * math.sin(lmda_rad), r_K])
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


# def rng_el(r):
#     """Get range and elevation from SEZ vector"""
#     rng = np.linalg.norm(r, axis=1)
#     el = np.arcsin(r[2] / rng)
#     el *= RAD2DEG
#     return el, rng


# def razel(r):
#     """Get range, azimuth, and elevation from SEZ vector"""
#     rng = np.linalg.norm(r, axis=0)
#     el = np.arcsin(r[2] / rng) * RAD2DEG
#     tmp = np.arctan2(r[0], r[1])
#     az = (tmp + np.pi * 0.5) * RAD2DEG
#     idx = np.all([r[0] < 0,r[1] < 0], axis=0)
#     az[idx] %= 360
#     return rng, az, el


# def azimuth(s, e):
#     """Compute azimuth from topocentric horizon coordinates SEZ
#     Args:
#         s : south vector from SEZ coordinate
#         e : east vector from SEZ coordinate
#     Output:
#         azimuth angle in radians with 0 at north.
#     """
#     out = np.arctan2(s, e) + math.pi * 0.5
#     if s < 0 and e < 0:
#         out = out % (2 * math.pi)
#     return out


# def elevation(z, rhomag):
#     """Compute elevation angle from topocentric horizon coordinates SEZ
#     Args:
#         z : Z vector from SEZ coordinate
#         rhomag : magnitude of SEZ coordinate vector
#     Output:
#         elevation angle in radians with 0 at horizon, pi/2 straight up
#     """
#     return np.arcsin(z / rhomag)


# def site_sat_rotations(rsiteECEF, rsatECEF):
#     rho = rsatECEF - np.array([[rsiteECEF[0]],[rsiteECEF[1]],[rsiteECEF[2]]], dtype=np.float64)
#     return rho