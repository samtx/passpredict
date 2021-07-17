# cython: boundscheck=False, wraparound=False

cimport cython

import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_f64_t

cdef extern from "passpredict.h":
    void c_sez2rngel(double *r_view, double *rng_view, double *el_view, int n)


def compute_range_and_elevation(np.ndarray[DTYPE_f64_t, ndim=2] r):
    """
    Compute range and elevation vectors from SEZ position vector
    """
    cdef int n = r.shape[0]
    cdef double[::1] r_view = r.flatten()
    cdef np.ndarray[DTYPE_f64_t, ndim=1] rng
    cdef np.ndarray[DTYPE_f64_t, ndim=1] el
    cdef double[::1] el_view
    cdef double[::1] rng_view
    rng_view = np.zeros(n, dtype=np.float64)
    el_view = np.zeros(n, dtype=np.float64)
    
    c_sez2rngel(&r_view[0], &rng_view[0], &el_view[0], n)
    rng = np.asarray(rng_view)
    el = np.asarray(el_view)
    return rng, el


      

