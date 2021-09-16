# cython: boundscheck=False, wraparound=False

cimport cython

cimport numpy as np
import numpy as np

ctypedef np.float64_t DTYPE_f64_t

cdef extern from "passpredict.h":
    void c_mod2ecef(double *jd, double *r, int n)
    void c_teme2ecef(double *jd, double *r, int n)
    void c_ecef2sez(double *r, double phi, double lmda, int n)


def mod2ecef(np.ndarray[np.float64_t, ndim=1] jd, np.ndarray[np.float64_t, ndim=2] r):
    """
    Rotate position vector r from MOD -> ECEF
    Using 1980 Nutation theory with GMST

    Inputs:
        jd, double[n]: Julian date floats
        r,  double[n, 3]: position vectors

    Outputs:
        r, double[n, 3]: rotated position vector rotated

    Reference:
        Uses IAU SOFA iauNut80 and iauGMST82 function to compute rotation matrix
    """
    cdef int n = jd.shape[0]
    cdef double[::1] jd_view = jd
    cdef double[::1] r_view = r.flatten()
    
    c_mod2ecef(&jd_view[0], &r_view[0], n)
    r_view_array = np.asarray(r_view)
    r = r_view_array.reshape((n, 3))
    return r


def teme2ecef(np.ndarray[np.float64_t, ndim=1] jd, np.ndarray[np.float64_t, ndim=2] r):
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
    cdef double[::1] r_view = r.flatten(order='C')
    
    c_teme2ecef(&jd_view[0], &r_view[0], n)
    r_view_array = np.asarray(r_view, dtype=np.float64)
    r = r_view_array.reshape((n, 3), order='C')
    return r


def ecef2sez(np.ndarray[np.float64_t, ndim=2] r, double phi, double lmda):
    """
    Rotate position vector r from ECEF -> SEZ
    
    Inputs:
        r,  double[n, 3]: position vectors
        phi, double: latitude in degrees
        lmda, double: longitude in degrees

    Outputs:
        r, double[n, 3]: rotated position vector rotated in-place

    """
    cdef int n = r.shape[0]
    cdef double[::1] r_view = r.flatten()
    
    c_ecef2sez(&r_view[0], phi, lmda, n)
    r_view_array = np.asarray(r_view)
    r = r_view_array.reshape((n, 3))
    return r



