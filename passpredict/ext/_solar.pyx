# cython: boundscheck=False, wraparound=False

cimport cython

import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_f64_t

cdef extern from "passpredict.h":
    void c_sun_pos(double *jd, double *r, int n)
    void c_sun_pos_ecef(double *jd, double *r, int n)
    void c_sun_sat_illumination_distance(double *rsat, double *rsun, double *illum_dist, int n)


def sun_pos(np.ndarray[DTYPE_f64_t, ndim=1] jd):
    """
    Get sun position vector in MOD frame
    """
    cdef int n = jd.shape[0]
    cdef double[::1] jd_view = jd
    cdef double[::1] r_view
    r_view = np.zeros(3*n, dtype=np.float64)
    
    c_sun_pos(&jd_view[0], &r_view[0], n)
    r_view_array = np.asarray(r_view)
    r = r_view_array.reshape((n, 3))
    return r


def sun_pos_ecef(np.ndarray[DTYPE_f64_t, ndim=1] jd):
    """
    Rotate position vector r from TEME -> ECEF (ITRF)
    Using 1976 Precession, 1980 Nutation theory, and 1982 GMST

    No polar motion, UT1, or EOP corrections

    Inputs:
        jd, double[n]: Julian date floats
        r,  double[n, 3]: position vectors

    Outputs:
        r, double[n, 3]: rotated position vector rotated in-place

    Reference:
        Uses IAU SOFA iauNut80 and iauGMST82 function to compute rotation matrix
    """
    cdef int n = jd.shape[0]
    cdef double[::1] jd_view = jd
    #cdef np.ndarray[DTYPE_f64_t, ndim=2] r
    cdef double[::1] r_view
    r_view = np.zeros(3*n, dtype=np.float64)
    
    c_sun_pos_ecef(&jd_view[0], &r_view[0], n)
    r = np.asarray(r_view).reshape((n, 3))
    return r
    
      
def sun_sat_illumination_distance(np.ndarray[DTYPE_f64_t, ndim=2] rsat, np.ndarray[DTYPE_f64_t, ndim=2] rsun):
    """
    Find illumination distance of satellite
    """
    cdef int n = rsat.shape[0]
    cdef double[::1] rsat_view = rsat.flatten()
    cdef double[::1] rsun_view = rsun.flatten()
    cdef np.ndarray[np.float64_t, ndim=1] illum_dist_array
    cdef double[::1] illum_dist_view = np.zeros(n, dtype=np.float64)

    c_sun_sat_illumination_distance(&rsat_view[0], &rsun_view[0], &illum_dist_view[0], n)
    illum_dist_array = np.asarray(illum_dist_view)
    return illum_dist_array

      

