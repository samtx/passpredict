# distutils: language = c

from libc.math cimport sin, cos, pi, asin, sqrt
import numpy as np
cimport numpy as np
cimport cython

from .constants import R_EARTH

from ._rotations import mod2ecef, mod2ecef_mjd


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void sun_pos_mod(double jd, double[::1] rmod):
    """
    Compute the Sun position vector.
    Return position vector in MOD coordinate frame
    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    """
    cdef double t_ut1, lmda_Msun, t_tdb, M_sun, lmda_eclp, lmda_eclp_rad, r_sun_mag, eps
    cdef double DEG2RAD = pi / 180.0   # degrees to radians
    cdef double AU_KM = 149597870.700  # AU to km
    r = np.empty(3, dtype=np.double)
    cdef double[::1] r_view = r

    jdld = np.longdouble(jd)
    t_ut1 = (jdld - 2451545.0) / 36525
    t_tdb = t_ut1
    # lmda_Msun = (280.4606184 + 36000.77005361 * t_tdb) % 360
    lmda_Msun = (280.460 + 36000.771*t_tdb)
    M_sun = (357.5291092 + 35999.05034*t_tdb) * DEG2RAD
    # M_sun = (357.5277233 + 35999.05034 * t_tdb) % 360
    lmda_eclp = lmda_Msun + 1.914666471*sin(M_sun) + 0.019994643*sin(2*M_sun)
    r_sun_mag = 1.000140612 - 0.016708617*cos(M_sun) - 0.000139589*cos(2*M_sun)
    eps = (23.439291 - 0.0130042*t_tdb) * DEG2RAD
    lmda_eclp_rad = lmda_eclp * DEG2RAD
    sinlmda = sin(lmda_eclp_rad)
    rmod[0] = r_sun_mag * cos(lmda_eclp_rad) * AU_KM
    rmod[1] = r_sun_mag * cos(eps) * sinlmda * AU_KM
    rmod[2] = r_sun_mag * sin(eps) * sinlmda * AU_KM


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sun_pos(double jd):
    """
    Compute the Sun position vector.
    Return position vector in ECEF coordinate frame
    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    """
    rmod = np.empty(3, dtype=np.double)
    recef = np.empty(3, dtype=np.double)
    cdef double[::1] rmod_view = rmod
    cdef double[::1] recef_view = recef

    sun_pos_mod(jd, rmod_view)
    mod2ecef(jd, rmod_view, recef_view)
    return recef


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sun_pos_mjd(double mjd):
    """
    Compute the Sun position vector.
    Return position vector in ECEF coordinate frame
    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    """
    rmod = np.empty(3, dtype=np.double)
    recef = np.empty(3, dtype=np.double)
    cdef double[::1] rmod_view = rmod
    cdef double[::1] recef_view = recef

    sun_pos_mod_mjd(mjd, rmod_view)
    mod2ecef_mjd(mjd, rmod_view, recef_view)
    return recef


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void sun_pos_mod_mjd(double mjd, double[::1] rmod):
    """
    Compute the Sun position vector.
    Return position vector in MOD coordinate frame
    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    """
    cdef double t_ut1, lmda_Msun, t_tdb, M_sun, lmda_eclp, lmda_eclp_rad, r_sun_mag, eps
    cdef double DEG2RAD = pi / 180.0   # degrees to radians
    cdef double AU_KM = 149597870.700  # AU to km
    r = np.empty(3, dtype=np.double)
    cdef double[::1] r_view = r

    # Find julian centuries since J2000 epoch
    # J2000 JD = 2451545.0
    # MJD JD = 2400000.5
    # mjd = jd - 2400000.5
    # jd = mjd + 2400000.5
    # --> jd - 2451545.0 --> (mjd + 2400000.5) - 2451545.0 --> mjd - 51544.5
    t_ut1 = (mjd - 51544.5) / 36525
    t_tdb = t_ut1
    # lmda_Msun = (280.4606184 + 36000.77005361 * t_tdb) % 360
    lmda_Msun = (280.460 + 36000.771*t_tdb)
    M_sun = (357.5291092 + 35999.05034*t_tdb) * DEG2RAD
    # M_sun = (357.5277233 + 35999.05034 * t_tdb) % 360
    lmda_eclp = lmda_Msun + 1.914666471*sin(M_sun) + 0.019994643*sin(2*M_sun)
    r_sun_mag = 1.000140612 - 0.016708617*cos(M_sun) - 0.000139589*cos(2*M_sun)
    eps = (23.439291 - 0.0130042*t_tdb) * DEG2RAD
    lmda_eclp_rad = lmda_eclp * DEG2RAD
    sinlmda = sin(lmda_eclp_rad)
    rmod[0] = r_sun_mag * cos(lmda_eclp_rad) * AU_KM
    rmod[1] = r_sun_mag * cos(eps) * sinlmda * AU_KM
    rmod[2] = r_sun_mag * sin(eps) * sinlmda * AU_KM

cpdef double sun_sat_angle(double[::1] rsat, double[::1] rsun):
    """Compute the sun-satellite angle
    Args:
        rsat : satellite position vector in ECI coordinates
        rsun : sun position vector in ECI coordinates
    Output:
        angle in radians between the two vectors
    References:
        Vallado, p. 912, Alg. 74
    """
    cdef double sun_cross_sat[3]
    cdef double[::1] sun_cross_sat_v = sun_cross_sat
    cdef double numer, denom, sinzeta, zeta

    # cross product
    # cprod[0] = a[1]*b[2] - a[2]*b[1]
    # cprod[1] = a[2]*b[0] - a[0]*b[2]
    # cprod[2] = a[0]*b[1] - a[1]*b[0]
    cross(rsun, rsat, sun_cross_sat_v)
    numer = norm(sun_cross_sat_v)
    denom = norm(rsun) * norm(rsat)
    sinzeta = numer / denom
    if (sinzeta > 1) and (1.0000003 > sinzeta):
        sinzeta = 1
    zeta = asin(sinzeta)
    return zeta


cpdef double sun_sat_orthogonal_distance(double[::1] rsat, double zeta):
    """
    Args:
        rsat : satellite position vector in ECI coordinates
        zeta : angle in radians between the satellite and sun vectors
    Output:
        distance from satellite to center of Earth orthogonal to sun vector
    """
    return norm(rsat) * cos(zeta - pi*0.5)


cpdef double sat_illumination_distance(double[::1] rsat, double[::1] rsun):
    cdef double zeta, dist
    zeta = sun_sat_angle(rsat, rsun)
    dist = sun_sat_orthogonal_distance(rsat, zeta)
    return dist


def is_sat_illuminated(double[::1] rsat, double[::1] rsun):
    dist = sat_illumination_distance(rsat, rsun)
    return dist > R_EARTH


cdef cross(double[::1] a, double[::1] b, double[::1] axb):
    """
    Cross product between two vectors
    """
    axb[0] = a[1]*b[2] - a[2]*b[1]
    axb[1] = a[2]*b[0] - a[0]*b[2]
    axb[2] = a[0]*b[1] - a[1]*b[0]
    return axb


cdef double norm(double[::1] a):
    """
    L2 norm of vector
    """
    cdef double n
    n = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    return n
