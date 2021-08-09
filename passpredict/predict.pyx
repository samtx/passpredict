# cython: boundscheck=False, wraparound=False
# cython: language_level=3

# Reference: https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html

cimport cython
from libcpp.list cimport list as stdlist
from libcpp.memory cimport shared_ptr
from .predict cimport (
    Point as Point_cpp,
    PassType as PassType_cpp,
    Overpass as Overpass_cpp,
    Orbit as Orbit_cpp,
    Location as Location_cpp,
    Satellite as Satellite_cpp,
    Observer as Observer_cpp,
    elsetrec
)
from .schemas import (
    Overpass as Overpass_py,
    Point as Point_py,
    PassType as PassType_py,
)


cdef extern from "ext/passpredict.h" namespace "passpredict":
    cdef Observer_cpp MakeObserver(Location_cpp, Satellite_cpp)

cdef class Point:
    cdef Point_cpp _point

    def __cinit__(self, double az, double el, double rng, double jd):
        self._point = Point_cpp()
        self._point.az = az
        self._point.el = el
        self._point.rng = rng
        self._point.jd = jd

    @property
    def jd(self):
        return self._point.jd

    @property
    def az(self):
        return self._point.az

    @property
    def el(self):
        return self._point.el

    @property
    def rng(self):
        return self._point.rng

cdef class Overpass:
    cdef Overpass_cpp _overpass

    def __cinit__(self):
        pass

cdef cppclass Observer:
    Observer_cpp* _observer

    def __cint__(self, Location_cpp location, Satellite_cpp satellite):
        self._observer = new Observer_cpp(location, satellite)

    def __dealloc__(self):
        del self._observer

    def get_overpasses(double t0, double tmax):
        overpasses = self._observer.GetOverpasses(t0, tmax)
        return overpasses


cdef predict(location_py, orbit_py, t0, tmax, min_el):
    """
    Create Observer cpp object from location and satellite and find predicted overpasses
    """
    cdef stdlist[shared_ptr[Overpass_cpp]] overpasses
    cdef int n, i
    cdef Location_cpp location
    location.lat_ = location_py.lat
    location.lon_ = location_py.lon
    location.h_ = location_py.h

    cdef Orbit_cpp orbit
    cdef elsetrec satrec
    satrec.satnum = orbit_py.satnum
    satrec.jdsatepoch = orbit_py.jdsatepoch
    satrec.jdsatepochF = orbit_py.jdsatepochF
    satrec.bstar = orbit_py.bstar
    satrec.inclo = orbit_py.inclo
    satrec.nodeo = orbit_py.nodeo
    satrec.ecco = orbit_py.ecco
    satrec.argpo = orbit_py.argpo
    satrec.mo = orbit_py.mo
    satrec.no_kozai = orbit_py.no_kozai
    orbit.satrec_ = satrec

    cdef Satellite_cpp satellite
    satellite.orbit_ = orbit

    observer_ptr = new Observer_cpp(location, satellite)
    try:
        overpasses = observer_ptr.GetOverpasses(t0, tmax)
    finally:
        del observer_ptr

    n = overpasses.size()

    # preallocate python list
    overpasses_py = [None]*n

    # convert cpp overpass object to python object
    for overpass in overpasses:
        o = overpass.get()
        overpass_py = Overpass_py(
            start_pt=Point_py(o.aos.)
        )

    return overpasses_py

# Call the predict() c++ function