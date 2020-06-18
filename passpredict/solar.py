# solar.py
import math

import numpy as np
from numpy.linalg import norm

from .constants import DEG2RAD, RAD2DEG, AU_KM, R_EARTH


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
    crossproduct = np.cross(rsun, rsat, axisa=0, axisb=0).T
    sinzeta = norm(crossproduct, axis=0) / (norm(rsun, axis=0) * norm(rsat, axis=0))
    return np.arcsin(sinzeta)


def sun_sat_orthogonal_distance(rsat, zeta):
    """
    Args:
        rsat : satellite position vector in ECI coordinates
        zeta : angle in radians between the satellite and sun vectors
    Output:
        distance from satellite to center of Earth orthogonal to sun vector
    """
    return norm(rsat, axis=0) * np.cos(zeta - math.pi * 0.5)


def is_sat_illuminated(rsat, rsun):
    zeta = sun_sat_angle(rsat, rsun)
    dist = sun_sat_orthogonal_distance(rsat, zeta)
    return dist > R_EARTH