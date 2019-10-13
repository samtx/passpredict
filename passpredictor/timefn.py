import datetime
import numpy as np
from .constants import DAY_S, J2000


def utc2tt(UTC, deltaAT=37., deltaUT1=0.):
    """Compute terrestial time from UTC

    deltaAT is posted annualy in the Astronomical Almanac. As of 2019 the offset
    is 37 seconds. deltaUT1 is posted daily by the US Navy

    Returns:
        TT : dynamic terrestial time in Julian Centuries

    References:
        Vallado, Section 3.5.4, Alg. 16
    """
    assert abs(deltaUT1) < 1
    UT1 = UTC + deltaUT1 / DAY_S
    # atomic time
    TAI = UTC + deltaAT / DAY_S
    # terrestial time in Julian Days
    JDtt = TAI + 32.184 / DAY_S
    # Convert JD to Julian Centuries
    TT = (JDtt - J2000) / 36525
    return TT


def julian_date(yr, mo=None, dy=None, hr=None, mn=None, sec=None):
    """Compute Julian Date from datetime or equivalent elements

    Notes
        T0 = 2451545.0

    References:
        Vallado, Algorithm 14, p.183
    """
    if isinstance(yr, datetime.datetime) or isinstance(yr, np.datetime64):
        if isinstance(yr, np.datetime64):
            dt = yr.astype(datetime.datetime)
        else:
            dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
        sec += dt.microsecond * (10**-6)
    jd1 = 367 * yr
    jd2 = 7 * (yr + (mo + 9) // 12) // 4
    jd3 = (275 * mo) // 9
    jd4 = dy
    jd5 = 1721013.5
    jd6 = ((sec / 60 + mn) / 60 + hr) / 24
    jd = jd1 - jd2 + jd3 + jd4 + jd5 + jd6
    # print([jd1, jd2, jd3, jd4, jd5, jd6])
    return jd


def julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day int."""
    janfeb = month < 3
    return (day
            + 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            - 32075)


def julian_date2(yr, mo=1, dy=1, hr=0, mn=0, sec=0.0):
    """Given a proleptic Gregorian calendar date, return a Julian date float."""
    if isinstance(yr, datetime.datetime) or isinstance(yr, np.datetime64):
        if isinstance(yr, np.datetime64):
            dt = yr.astype(datetime.datetime)
        else:
            dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
        sec += dt.microsecond * (10**-6)

    return julian_day(yr, mo, dy) - 0.5 + (
        sec + mn * 60.0 + hr * 3600.0) / DAY_S
