# cython: boundscheck=False, wraparound=False
# cython: language_level=3
# distutils: language = c

from libc.math cimport sin, cos, sqrt, atan, atan2, asin, pi, pow, fmod

import numpy as np
cimport numpy as np


cdef extern from "sofa.h":
    cdef double iauGmst82(double dj1, double dj2)
    cdef double iauEqeq94(double date1, double date2)
    cdef double iauObl80(double date1, double date2)
    cdef void iauNut80(double date1, double date2, double *dpsi, double *deps)
    cdef void iauNumat(double epsa, double dpsi, double deps, double rmatn[3][3])
    cdef int iauJd2cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd)
    cdef int iauDat(int iy, int im, int id, double fd, double *deltat)
    cdef int iauUtctai(double utc1, double utc2, double *tai1, double *tai2)
    cdef int iauTaitt(double tai1, double tai2, double *tt1, double *tt2)
    cdef int iauD2dtf(const char *scale, int ndp, double d1, double d2, int *iy, int *im, int *id, int ihmsf[4])
    cdef void iauIr(double r[3][3])
    cdef void iauRz(double psi, double r[3][3])
    cdef void iauTr(double r[3][3], double rt[3][3])
    cdef void iauRxr(double a[3][3], double b[3][3], double atb[3][3])
    cdef void iauRxp(double r[3][3], double p[3], double rp[3])


cpdef void ecef_to_rhosez(
    double location_lat_rad,
    double location_lon_rad,
    double[::1] location_ecef,
    double[::1] satellite_ecef,
    double[::1] rho_sez,
):
    """
    Get slant vector to satellite relative for observing location
    """
    cdef double rx, ry, rz
    cdef double sin_location_lat, cos_location_lat
    cdef double sin_location_lon, cos_location_lon

    rx = satellite_ecef[0] - location_ecef[0]
    ry = satellite_ecef[1] - location_ecef[1]
    rz = satellite_ecef[2] - location_ecef[2]

    sin_location_lat = sin(location_lat_rad)
    cos_location_lat = cos(location_lat_rad)
    sin_location_lon = sin(location_lon_rad)
    cos_location_lon = cos(location_lon_rad)

    rho_sez[0] = (sin_location_lat * cos_location_lon * rx) + (sin_location_lat * sin_location_lon * ry) - (cos_location_lat * rz)
    rho_sez[1] = -sin_location_lon * rx + cos_location_lon * ry
    rho_sez[2] = (cos_location_lat * cos_location_lon * rx) + (cos_location_lat * sin_location_lon * ry) + (sin_location_lat * rz)


def range_at(
    double location_lat_rad,
    double location_lon_rad,
    double[::1] location_ecef,
    double[::1] satellite_ecef
):
    """
    Get slant range of satellite relative for observing location [km]
    """
    cdef double rho[3]
    cdef double[::1] rho_v = rho
    cdef double range_
    ecef_to_rhosez(location_lat_rad, location_lon_rad, location_ecef, satellite_ecef, rho_v)
    range_ = sqrt(rho_v[0]*rho_v[0] + rho_v[1]*rho_v[1] + rho_v[2]*rho_v[2])
    return range_


def elevation_at(
    double location_lat_rad,
    double location_lon_rad,
    double[::1] location_ecef,
    double[::1] satellite_ecef
):
    """
    Get elevation of satellite relative for observing location [deg]
    """
    cdef double rho[3]
    cdef double[::1] rho_v = rho
    cdef double range_, el, el_deg
    ecef_to_rhosez(location_lat_rad, location_lon_rad, location_ecef, satellite_ecef, rho_v)
    range_ = sqrt(rho_v[0]*rho_v[0] + rho_v[1]*rho_v[1] + rho_v[2]*rho_v[2])
    el = asin(rho_v[2] / range_)
    # convert radians to degrees
    el_deg = el * 180.0 / pi
    return el_deg


def elevation_at_rad(
    double cos_lat_x_cos_lon,
    double cos_lat_x_sin_lon,
    double sin_lat,
    double[::1] location_ecef,
    double[::1] satellite_ecef
):
    """
    Get elevation of satellite relative for observing location [rad]
    """
    cdef double rho_z
    cdef double range_, el
    cdef double rx, ry, rz
    rx = satellite_ecef[0] - location_ecef[0]
    ry = satellite_ecef[1] - location_ecef[1]
    rz = satellite_ecef[2] - location_ecef[2]
    rho_z = (cos_lat_x_cos_lon * rx) + (cos_lat_x_sin_lon * ry) + (sin_lat * rz)
    range_ = sqrt(rx*rx + ry*ry + rz*rz)
    el = asin(rho_z / range_)
    return el


def razel(
    double location_lat_rad,
    double location_lon_rad,
    double[::1] location_ecef,
    double[::1] satellite_ecef,
):
    """
    Get range, azimuth, and elevation of satellite relative for observing location
    """
    cdef double rho[3]
    cdef double range_, el, az, el_deg, az_deg
    ecef_to_rhosez(location_lat_rad, location_lon_rad, location_ecef, satellite_ecef, rho)
    range_ = sqrt(rho[0]*rho[0] + rho[1]*rho[1] + rho[2]*rho[2])
    el = asin(rho[2] / range_)
    az = atan2(-rho[1], rho[0]) + pi
    # convert radians to degrees
    el_deg = el * 180.0 / pi
    az_deg = az * 180.0 / pi
    return (range_, az_deg, el_deg)


def ecef_to_llh(double[:] recef):
    """
    Convert ECEF coordinates to latitude, longitude, and altitude
    Uses WGS84 constants
    Based on orbit_predictor.coordinate_systems.ecef_to_llh
    """
    # WGS-84 ellipsoid parameters */
    cdef double a = 6378.1370
    cdef double b = 6356.752314
    cdef double p, thet, esq, epsq, lat, lon, h, n

    p = sqrt(recef[0]*recef[0] + recef[1]*recef[1])
    thet = atan(recef[2] * a / (p * b))
    esq = 1.0 - (b / a)*(b / a)
    epsq = (a / b)*(a / b) - 1.0

    lat = atan((recef[2] + epsq * b * pow(sin(thet), 3)) / (p - esq * a * pow(cos(thet), 3)))
    lon = atan2(recef[1], recef[0])
    n = a*a / sqrt(a*a*cos(lat)*cos(lat) + b*b*sin(lat)*sin(lat))
    h = p / cos(lat) - n

    lat = lat * 180.0 / pi  # convert from radians to degrees
    lon = lon * 180.0 / pi
    return lat, lon, h


cpdef double jd2tt(double jd):
    """
    Convert julian date to terrestial time. Don't apply corrections for UT1
    """
    cdef double tt1
    cdef double tt2
    cdef double tt
    # Find terrestial time, ignore delta_UT1
    jd2tt2(jd, &tt1, &tt2)
    tt = tt1 + tt2
    return tt


cdef void jd2tt2(double jd, double* tt1, double* tt2):
    """
    Convert julian date to terrestial time. Don't apply corrections for UT1
    """
    cdef double tai1, tai2
    cdef int err
    # Find terrestial time, ignore delta_UT1
    err = iauUtctai(jd, 0.0, &tai1, &tai2)
    err = iauTaitt(tai1, tai2, tt1, tt2)
    return


cdef void mjd2tt2(double mjd, double* tt1, double* tt2):
    """
    Convert modified julian date to terrestial time. Don't apply corrections for UT1
    """
    cdef double tai1, tai2
    cdef int err
    # Find terrestial time, ignore delta_UT1
    err = iauUtctai(2400000.5, mjd, &tai1, &tai2)
    err = iauTaitt(tai1, tai2, tt1, tt2)
    return


cpdef mod2ecef(double jd, double[::1] rmod, double[::1] recef):
    """
    Convert MOD to ECEF coordinates

    N = nutation rotation matrix
    G = z-rotation matrix by gast

    r_mod = [NG]r_ecef
    --> r_ecef = [NG]^T r_mod

    """
    cdef double dp80, de80, epsa, tt1, tt2, ee, gast
    cdef double N[3][3]
    cdef double G[3][3]
    cdef double NG[3][3]
    cdef double NGT[3][3]
    cdef double twopi = 2*pi

    # get terrestial time
    jd2tt2(jd, &tt1, &tt2)
    # get nutation values
    iauNut80(tt1, tt2, &dp80, &de80)
    # mean obliquity
    epsa = iauObl80(tt1, tt2)
    # build nutation rotation matrix
    iauNumat(epsa, dp80, de80, N)
    # equation of equinoxes
    ee = iauEqeq94(tt1, tt2)
    # greenwich apparent sidereal time
    gast = iauGmst82(jd, 0.0) + ee
    # normalize gast into 0 <= gast < 2pi
    gast = fmod(gast, twopi)
    if gast < 0:
        gast += twopi
    iauIr(G)  # initialize G matrix with identity
    iauRz(gast, G)   # rotate on the Z axis by gast radians
    """
    Rotate on Z axis
    (  + cos(psi)   + sin(psi)     0  )
    (                                 )
    (  - sin(psi)   + cos(psi)     0  )
    (                                 )
    (       0            0         1  )
    """

    # Create rotation matrix, multiply NG, then transpose
    iauRxr(N, G, NG)
    # iauTr(NG, NG)
    iauRxp(NG, &rmod[0], &recef[0])


cpdef mod2ecef_mjd(double mjd, double[::1] rmod, double[::1] recef):
    """
    Convert MOD to ECEF coordinates

    N = nutation rotation matrix
    G = z-rotation matrix by gast

    r_mod = [NG]r_ecef
    --> r_ecef = [NG]^T r_mod

    """
    cdef double dp80, de80, epsa, tt1, tt2, ee, gast
    cdef double N[3][3]
    cdef double G[3][3]
    cdef double NG[3][3]
    cdef double NGT[3][3]
    cdef double twopi = 2*pi

    # get terrestial time
    mjd2tt2(mjd, &tt1, &tt2)
    # get nutation values
    iauNut80(tt1, tt2, &dp80, &de80)
    # mean obliquity
    epsa = iauObl80(tt1, tt2)
    # build nutation rotation matrix
    iauNumat(epsa, dp80, de80, N)
    # equation of equinoxes
    ee = iauEqeq94(tt1, tt2)
    # greenwich apparent sidereal time
    gast = iauGmst82(2400000.5, mjd) + ee
    # normalize gast into 0 <= gast < 2pi
    gast = fmod(gast, twopi)
    if gast < 0:
        gast += twopi
    iauIr(G)  # initialize G matrix with identity
    iauRz(gast, G)   # rotate on the Z axis by gast radians
    """
    Rotate on Z axis
    (  + cos(psi)   + sin(psi)     0  )
    (                                 )
    (  - sin(psi)   + cos(psi)     0  )
    (                                 )
    (       0            0         1  )
    """

    # Create rotation matrix, multiply NG, then transpose
    iauRxr(N, G, NG)
    # iauTr(NG, NG)
    iauRxp(NG, &rmod[0], &recef[0])


cpdef double gmst82(double mjd):
    """
    Get Greenwich Mean Sidereal time, 1982 model
    Input with modified julian date
    """
    cdef double _gmst
    _gmst = iauGmst82(2400000.5, mjd)
    return _gmst



def teme2ecef(double mjd, double[::1] rteme, double[::1] recef):
    """
    Convert a TEME position vector to ECEF.
    Inputs modified julian date
    """
    cdef double gmst, singmst, cosgmst
    gmst = gmst82(mjd)
    singmst = sin(gmst)
    cosgmst = cos(gmst)
    recef[0] = rteme[0]*cosgmst + rteme[1]*singmst
    recef[1] = rteme[0]*(-singmst) + rteme[1]*cosgmst
    recef[2] = rteme[2]
