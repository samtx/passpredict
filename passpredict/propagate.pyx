
from propagate cimport pkepler as pkepler_cpp

def pkepler(double[:] r0, double[:] v0, double dtsec, double ndot = 0, double nddot = 0):
    """
    pkepler algorithm from Vallado
    """
    cdef double[3] r = [0, 0, 0]
    cdef double[3] v = [0, 0, 0]
    pkepler_cpp(&r0[0], &v0[0], dtsec, ndot, nddot, &r[0], &v[0])
    return (r, v)

def sgp4():
    pass
