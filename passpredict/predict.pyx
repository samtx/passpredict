# cython: boundscheck=False, wraparound=False

cimport cython

import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_f64_t

cdef extern from "ext/passpredict.hpp":
    void c_sun_pos(double *jd, double *r, int n)
    void c_sun_pos_ecef(double *jd, double *r, int n)
    void c_sun_sat_illumination_distance(double *rsat, double *rsun, double *illum_dist, int n)

# Call the predict() c++ function