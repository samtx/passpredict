# cython: boundscheck=False, wraparound=False
# cython: language_level=3

# Reference: https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html

cimport cython
from libcpp.list cimport list as stdlist
from libcpp.memory cimport shared_ptr

import numpy as np
cimport numpy as np

from collections import namedtuple

from passpredict.predictcpp cimport (
    Location as Location_cpp,
)
from passpredict import OMM


# cdef extern from "ext/passpredict.h" namespace "passpredict":
#     # cdef Observer_cpp MakeObserver(Location_cpp, Satellite_cpp)
#     cdef double ComputeElevationAngle(double, double, Location_cpp, Satellite_cpp)


cdef class Location:
    """
    Python wrapper for passpredict::Location cpp class
    """
    cdef Location_cpp* _location

    def __cinit__(self, double lat, double lon, double h, str name=""):
        self._location = new Location_cpp(lat, lon, h)
        self.name = name

    def __dealloc__(self):
        del self._location

    @property
    def lat(self):
        return self._location.lat_

    @property
    def lon(self):
        return self._location.lon_

    @property
    def h(self):
        return self._location.h_


cdef class Satellite:
    """
    Python wrapper for passpredict::Satellite cpp class
    """
    cdef Satellite_cpp* _satellite

    def __cinit__(self, ):
        self.brightness = None

    def __dealloc__(self):
        pass

    @classmethod
    def from_tle(cls, str tle1, str tle2, brightness = None):
        cdef Satellite satellite = Satellite.__new__(self)




# cdef class Point:
#     cdef Point_cpp _point

#     def __cinit__(self, double az, double el, double rng, double jd):
#         self._point = Point_cpp()
#         self._point.az = az
#         self._point.el = el
#         self._point.rng = rng
#         self._point.jd = jd

#     @property
#     def jd(self):
#         return self._point.jd

#     @property
#     def az(self):
#         return self._point.az

#     @property
#     def el(self):
#         return self._point.el

#     @property
#     def rng(self):
#         return self._point.rng


# # cdef class Overpass:
# #     cdef Overpass_cpp _overpass

# #     def __cinit__(self):
# #         pass


# ctypedef fused array:
#     double


# class to hold azimuth, elevation, and range results
AzElRng = namedtuple('AzElRng', 'az el rng')

# cdef class Observer:
#     cdef Observer_cpp* _observer

#     def __cint__(self, Location location, Satellite_cpp satellite):
#         self._observer = new Observer_cpp(location._location, satellite)

#     def __dealloc__(self):
#         del self._observer

#     def get_overpasses(self, double t0, double tmax):

#         overpasses = self._observer.GetOverpasses(t0, tmax)
#         # return overpasses

#     def compute_azelrng(self, t):
#         # return array

#         cdef int i, n

#         if hasattr(t, "__len__"):
#             # compute for single time step
#             self._observer.UpdateToJd(t)
#             result = AzElRng(
#                 self._observer.az_,
#                 self._observer.el_,
#                 self._observer.range_
#             )
#             return result

#         # compute for each time step
#         n = t.size()
#         cdef double[::1] az = np.zeros(n, dtype=np.double)
#         cdef double[::1] el = np.zeros(n, dtype=np.double)
#         cdef double[::1] rng = np.zeros(n, dtype=np.double)
#         for i in range(n):
#             self._observer.UpdateToJd(t[i])
#             az[i] = self._observer.az_
#             el[i] = self._observer.el_
#             rng[i] = self._observer.range_
#         result = AzElRng(az, el, rng)
#         return result



# cdef predict(location_py, orbit_py, t0, tmax, min_el):
#     """
#     Create Observer cpp object from location and satellite and find predicted overpasses
#     """
#     cdef stdlist[shared_ptr[Overpass_cpp]] overpasses
#     cdef int n, i
#     cdef Location_cpp location
#     location.lat_ = location_py.lat
#     location.lon_ = location_py.lon
#     location.h_ = location_py.h

#     cdef Orbit_cpp orbit
#     cdef elsetrec satrec
#     satrec.satnum = orbit_py.satnum
#     satrec.jdsatepoch = orbit_py.jdsatepoch
#     satrec.jdsatepochF = orbit_py.jdsatepochF
#     satrec.bstar = orbit_py.bstar
#     satrec.inclo = orbit_py.inclo
#     satrec.nodeo = orbit_py.nodeo
#     satrec.ecco = orbit_py.ecco
#     satrec.argpo = orbit_py.argpo
#     satrec.mo = orbit_py.mo
#     satrec.no_kozai = orbit_py.no_kozai
#     orbit.satrec_ = satrec

#     cdef Satellite_cpp satellite
#     satellite.orbit_ = orbit

#     observer_ptr = new Observer_cpp(location, satellite)
#     try:
#         overpasses = observer_ptr.GetOverpasses(t0, tmax)
#     finally:
#         del observer_ptr

#     n = overpasses.size()

#     # preallocate python list
#     overpasses_py = [None]*n

#     # convert cpp overpass object to python object
#     # for overpass in overpasses:
#     #     o = overpass.get()
#     #     overpass_py = Overpass_py(
#     #         start_pt=Point_py(o.aos.)
#     #     )

#     return overpasses_py

# Call the predict() c++ function

# def compute_elevation_angle(jd, location, satellite):
#     """
#     Compute the elevation angle from the horizon for an observer
#     """
#     if location.