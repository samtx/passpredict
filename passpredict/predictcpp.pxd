# cython: language_level=3

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

    cdef struct Omm:
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
        double mo    # deg
        double no_kozai

    cdef cppclass Orbit:
        Orbit()
        Orbit(Omm)
        elsetrec satrec_

    cdef cppclass Location:
        Location()
        Location(double, double, double)
        double lat_
        double lon_
        double h_
        double[3] recef_

    cdef cppclass Satellite:
        Satellite()
        Satellite(Orbit)
        Orbit orbit_