# cython: language_level=3

from libcpp.list cimport list as stdlist
from libcpp.memory cimport shared_ptr
from libcpp.string cimport string

cdef extern from "ext/observer.cpp":
    pass


cdef extern from "ext/passpredict.h":
    cdef struct elsetrec:
        char classification
        char satnum[9]
        int revnum
        int elnum
        int ephtype
        double jdsatepoch
        double jdsatepochF
        double bstar
        double inclo  # deg
        double nodeo  # deg
        double ecco
        double argpo  # deg
        double mo     # deg
        double no_kozai


cdef extern from "ext/passpredict.h" namespace "passpredict":
    # cdef cppclass Point:
    #     double az
    #     double el
    #     double rng
    #     double jd
    #     Point()

    # cdef enum PassType:
    #     visible = 0
    #     unlit = 1
    #     daylight = 2

    # cdef cppclass Overpass:
    #     PassType pass_type
    #     double brightness
    #     double altitude
    #     double duration
    #     double duration_vis
    #     Point aos
    #     Point los
    #     Point max
    #     Point aos_vis
    #     Point los_vis

    # cdef cppclass Orbit:
    #     Orbit()
    #     elsetrec satrec_

    cdef cppclass Location:
        Location()
        Location(double, double, double)
        double lat_
        double lon_
        double h_

    # cdef cppclass Satellite:
    #     Satellite()
    #     Orbit orbit_

    # cdef cppclass Observer:
    #     Observer(Location, Satellite)
    #     stdlist[shared_ptr[Overpass]] GetOverpasses(double t0, double tmax)
    #     void UpdateToJd(double)
    #     double az_, el_, range_
