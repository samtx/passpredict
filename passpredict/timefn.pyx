# cython: boundscheck=False, wraparound=False
# cython: language_level=3

cimport cython

cdef extern from "SGP4.h" namespace "SGP4Funcs":
    cdef void days2mdhms_SGP4(int year, double days, int& mon, int& day, int& hr, int& minute, double& sec)
    cdef void jday_SGP4(int year, int mon, int day, int hr, int minute, double sec, double& jd, double& jdFrac)

def epoch_to_jd(int epoch_year, double epoch_days):
    """
    Convert TLE epoch value to julian date
    """
    cdef int year, mon, day, hr, minute
    cdef double sec, jd, jdfr
    if epoch_year < 57:
        year = epoch_year + 2000
    else:
        year = epoch_year + 1900
    days2mdhms_SGP4(year, epoch_days, mon, day, hr, minute, sec)
    jday_SGP4(year, mon, day, hr, minute, sec, jd, jdfr)
    return (jd, jdfr)

def julian_date(int year, int mon, int day, int hr, int minute, double sec):
    """
    Convert date to julian date
    """
    cdef double jd
    jd = jday_SGP4(year, mon, day, hr, minute, sec, jd, jdfr)
    return jd





