# cython: boundscheck=False, wraparound=False
# cython: language_level=3

from libc.math cimport sin, cos, sqrt, atan2, asin, pi

import numpy as np
cimport numpy as np


def razel(
    double location_lat_rad,
    double location_lon_rad,
    double[:] location_ecef,
    double[:] satellite_ecef
):
    """
    Get range, azimuth, and elevation of satellite relative for observing location
    """
    cdef double rx, ry, rz
    cdef double sin_location_lat, cos_location_lat
    cdef double sin_location_lon, cos_location_lon
    cdef double top_s, top_e, top_z
    cdef double range_, el, az, el_deg, az_deg

    rx = satellite_ecef[0] - location_ecef[0]
    ry = satellite_ecef[1] - location_ecef[1]
    rz = satellite_ecef[2] - location_ecef[2]

    sin_location_lat = sin(location_lat_rad)
    cos_location_lat = cos(location_lat_rad)
    sin_location_lon = sin(location_lon_rad)
    cos_location_lon = cos(location_lon_rad)

    top_s = (sin_location_lat * cos_location_lon * rx) + (sin_location_lat * sin_location_lon * ry) - (cos_location_lat * rz)
    top_e = -sin_location_lon * rx + cos_location_lon * ry
    top_z = (cos_location_lat * cos_location_lon * rx) + (cos_location_lat * sin_location_lon * ry) + (sin_location_lat * rz)

    range_ = sqrt(top_s*top_s + top_e*top_e + top_z*top_z)
    el = asin(top_z / range_)
    az = atan2(-top_e, top_s) + pi

    # convert radians to degrees
    el_deg = el * 180.0 / pi
    az_deg = az * 180.0 / pi

    return (range_, az_deg, el_deg)