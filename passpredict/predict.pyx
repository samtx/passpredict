# cython: boundscheck=False, wraparound=False
# cython: language_level=3

# Reference: https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html

cimport cython
from libcpp.list cimport list as stdlist
from libcpp.memory cimport shared_ptr
from libc.string cimport strcpy

import numpy as np
cimport numpy as np

from collections import namedtuple

from .predictcpp cimport Location as Location_cpp
from .predictcpp cimport Omm as Omm_cpp
from .predictcpp cimport Orbit as Orbit_cpp
from .predictcpp cimport Satellite as Satellite_cpp
from passpredict.timefn import epoch_to_jd
from passpredict import OMM


# cdef extern from "ext/passpredict.h" namespace "passpredict":
#     # cdef Observer_cpp MakeObserver(Location_cpp, Satellite_cpp)
#     cdef double ComputeElevationAngle(double, double, Location_cpp, Satellite_cpp)


# class to hold azimuth, elevation, and range results
AzElRng = namedtuple('AzElRng', 'az el rng')


cdef class Location:
    """
    Python wrapper for passpredict::Location cpp class
    """
    cdef Location_cpp _location
    cdef str name

    def __cinit__(self, double lat, double lon, double h, name = None):
        self._location = Location_cpp(lat, lon, h)
        self.name = name

    @property
    def lat(self):
        return self._location.lat_

    @property
    def lon(self):
        return self._location.lon_

    @property
    def h(self):
        return self._location.h_

    @property
    def name(self):
        return self.name

    @property
    def recef(self):
        return self._location.recef_



cdef class Satellite:
    """
    Python wrapper for passpredict::Satellite cpp class
    """
    cdef Satellite_cpp _satellite
    cdef Orbit_cpp _orbit
    cdef Omm_cpp _omm

    def __cinit__(self, omm, brightness = None):
        # create cpp orbit object, then use that to create cpp satellite object
        self._omm.jdsatepoch = omm.jdsatepoch
        self._omm.jdsatepochF = omm.jdsatepochF
        self._omm.no_kozai = omm.no_kozai
        self._omm.ecco = omm.ecco
        self._omm.inclo = omm.inclo
        self._omm.nodeo = omm.nodeo
        self._omm.argpo = omm.argpo
        self._omm.mo = omm.mo
        self._omm.nddot = omm.nddot
        self._omm.bstar = omm.bstar
        self._omm.ndot = omm.ndot
        self._omm.elnum = omm.elnum
        self._omm.revnum = omm.revnum
        self._omm.classification = omm.classification
        self._omm.ephtype = omm.ephtype
        self._orbit = Orbit_cpp(self._omm)
        self._satellite = Satellite_cpp(self._orbit)
        self.brightness = None

    # def __dealloc__(self):
    #     del self._satellite
    #     del self._orbit





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

# cdef class Orbit:
#     cdef char classification
#     cdef char[9] satnum
#     cdef double jdsatepoch
#     cdef int revnum
#     cdef int elnum
#     cdef int ephtype
#     cdef double jdsatepochF
#     cdef double bstar
#     cdef double ecco
#     cdef double inclo
#     cdef double nodeo
#     cdef double argpo
#     cdef double mo
#     cdef double no_kozai

#     @staticmethod
#     cdef Orbit from_tle(unsigned char[:] tle1, unsigned char[:] tle2):
#         """
#         Convert TLE strings to Orbit object
#         """
#         satnum = tle1[2:7]
#         classification = tle1[8]
#         epoch_year = int(tle1[18:20])
#         epoch_days = float(tle1[20:32])
#         jdsatepoch, jdsatepochF = epoch_to_jd(epoch_year, epoch_days)
#         # ndot = float(tle1[34:44])
#         # nddot = float(tle1[45:52])
#         bstar = float(tle1[54:62])
#         ephtype = tle1[63]
#         elnum = int(tle1[65:69])
#         inclo = float(tle2[9:17])  # inclination
#         nodeo = float(tle2[18:26])  # right ascension of ascending node
#         ecco = float(tle2[27:34]) / 1e7  # eccentricity
#         argpo = float(tle2[35:43])
#         mo = float(tle2[44:52])    # mean anomaly
#         no_kozai = float(tle2[53:64])   # mean motion
#         revnum = int(tle2[64:69])
#         cdef Orbit orbit = cls(
#             satnum=satnum,
#             jdsatepoch=jdsatepoch,
#             jdsatepochF=jdsatepochF,
#             bstar=bstar,
#             inclo=inclo,
#             nodeo=nodeo,
#             ecco=ecco,
#             argpo=argpo,
#             mo=mo,
#             no_kozai=no_kozai,
#             revnum=revnum,
#             elnum=elnum,
#             classification=classification,
#             ephtype=ephtype
#         )
#         return orbit

#     @staticmethod
#     cdef Orbit from_omm(omm):
#         """
#         Convert OMM object to Orbit object
#         """
#         cdef Orbit orbit
#         strcpy(omm.satnum, orbit.satnum)
#         orbit.jdsatepoch = omm.jdsatepoch
#         orbit.jdsatepochF = omm.jdsatepochF
#         orbit.bstar = omm.bstar
#         orbit.inclo = omm.inclo
#         orbit.nodeo = omm.nodeo
#         orbit.ecco = omm.ecco
#         orbit.argpo = omm.argpo
#         orbit.mo = omm.mo
#         orbit.no_kozai = omm.no_kozai
#         orbit.revnum = omm.revnum
#         orbit.elnum = omm.elnum
#         strcpy(omm.classification, orbit.classification)
#         strcpy(omm.ephtype, orbit.ephtype)
#         return orbit
