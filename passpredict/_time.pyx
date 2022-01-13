# cython: boundscheck=False, wraparound=False
# cython: language_level=3
# distutils: language = c++

from libc.math cimport round as c_round, fmod, floor, fabs

import datetime

cdef extern from "SGP4.h" namespace "SGP4Funcs":
    cdef void days2mdhms_SGP4(int year, double days, int& mon, int& day, int& hr, int& minute, double& sec)
    cdef void jday_SGP4(int year, int mon, int day, int hr, int minute, double sec, double& jd, double& jdFrac)
    cdef void invjday_SGP4(double jd, double jdFrac, int& year, int& mon, int& day, int& hr, int& minute, double& sec)


tz_utc = datetime.timezone.utc


def epoch_to_jd(int epoch_year, double epoch_days):
    """
    Convert TLE epoch value to julian date
    """
    cdef int year, mon, day, hr, minute
    cdef double sec
    cdef double jd, jdfr
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
    cdef double jdfr
    jday_SGP4(year, mon, day, hr, minute, sec, jd, jdfr)
    return (jd, jdfr)


def jday2datetime(double jd):
    """
    Convert julian date float to python datetime object rounded to nearest second
    """
    cdef int year, mon, day, hr, minute, i_sec
    cdef double d_sec
    cdef double jdFrac = 0
    invjday_SGP4(jd, jdFrac, year, mon, day, hr, minute, d_sec)
    i_sec = int(c_round(d_sec))
    if i_sec == 60:
        return datetime.datetime(year, mon, day, hr, minute, 59, tzinfo=tz_utc) + datetime.timedelta(seconds=1)
    return datetime.datetime(year, mon, day, hr, minute, i_sec, tzinfo=tz_utc)


def jday2datetime_us(double jd):
    """
    Convert julian date float to python datetime object including microseconds
    """
    cdef int year, mon, day, hr, minute
    cdef double d_sec
    cdef double jdFrac = 0
    invjday_SGP4(jd, jdFrac, year, mon, day, hr, minute, d_sec)
    i_sec, i_us = divmod(d_sec, 1)
    return datetime.datetime(year, mon, day, hr, minute, int(i_sec), int(i_us*1e6), tzinfo=tz_utc)
