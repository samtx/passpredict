import datetime

import numpy as np
from astropy.time import Time

from .constants import DAY_S, J2000

tz_utc = datetime.timezone.utc


def jd2jc(jd1: float, jd2: float = 0.0) -> float:
    """
    Convert julian date to julian century
    """
    return ((jd1 - J2000) + jd2) / 36525.0


def utc2tt(UTC, deltaAT=37.0, deltaUT1=0.0):
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
    TT = jd2jc(JDtt)
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
        sec += dt.microsecond * (10 ** -6)
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
    return (
        day
        + 1461 * (year + 4800 - janfeb) // 4
        + 367 * (month - 2 + janfeb * 12) // 12
        - 3 * ((year + 4900 - janfeb) // 100) // 4
        - 32075
    )


def julian_date2(yr, mo=1, dy=1, hr=0, mn=0, sec=0.0):
    """Given a proleptic Gregorian calendar date, return a Julian date float."""
    if isinstance(yr, datetime.datetime) or isinstance(yr, np.datetime64):
        if isinstance(yr, np.datetime64):
            dt = yr.astype(datetime.datetime)
        else:
            dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
        sec += dt.microsecond * (10 ** -6)

    return julian_day(yr, mo, dy) - 0.5 + (sec + mn * 60.0 + hr * 3600.0) / DAY_S


def days2mdhms(year: int, days: float):
    """This procedure converts the day of the year, days, to the equivalent
    month, day, hour, minute and second.
    Args:
        year : int, year between 1900 - 2100
        days : float, julian day of the year between 0.0 - 366.0
    Outputs:
        mon : int, month, 1 .. 12
        day : int, day, 1 .. 28,29,30,31
        hr : int, hour, 0 .. 23
        min : int, minute 0 .. 59
        sec : float, second, 0.0 .. 59.999
    References:
        Vallado
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # non vectorized version...
    lmonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    dayofyr = int(days // 1.0)  # day of year
    # find month and day of month
    if (year % 4) == 0:
        # int array containing the number of days per month
        lmonth = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    i = 1
    inttemp = 0
    while dayofyr > inttemp + lmonth[i - 1] and i < 12:
        inttemp = inttemp + lmonth[i - 1]
        i += 1
    mon = i
    day = dayofyr - inttemp
    # find hours minutes and seconds
    temp = (days - dayofyr) * 24.0
    hr = int(temp // 1.0)
    temp = (temp - hr) * 60.0
    minute = int(temp // 1.0)
    sec = (temp - minute) * 60.0
    return mon, day, hr, minute, sec


def invjday(jd: float):
    """This procedure finds the year, month, day, hour, minute and second given the
    julian date. jd can be ut1, tdt, tdb, etc.
    Args:
        jd : float, julian date, days from 4713 BCE
    Outputs:
        year : int, year between 1900 - 2100
        mon : int, month between 1 - 12
        day : int, day between 1 - 31
        hr : int, hour between 0 - 23
        min : int, minute between 0 - 59
        sec : float, second between 0.0 - 59.999
    References:
        Vallado, 2007, 208, alg 22, ex 3-13
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # find year and days of the year
    temp = jd - 2415019.5
    tu = temp / 365.25  # julian centuries from 0 h jan 0, 1900
    year = 1900 + int(tu // 1.0)
    leapyrs = int(((year - 1901) * 0.25) // 1.0)

    # optional nudge by 8.64x10-7 sec to get even outputs
    days  = temp - ((year - 1900) * 365.0 + leapyrs) + 0.00000000001
    # check for case of beginning of a year
    if (days < 1.0):
        year    = year - 1
        leapyrs = int(((year - 1901) * 0.25) // 1.0)
        days    = temp - ((year - 1900) * 365.0 + leapyrs)
    # find remaing data
    mon, day, hr, minute, sec = days2mdhms(year, days)
    sec -= 0.00000086400
    return year, mon, day, hr, minute, sec


def jday2datetime(jdt: float) -> datetime.datetime:
    """Turn julian day into datetime object rounded to nearest second"""
    datetuple = invjday(jdt)
    yr, mo, date, hr, mn, sec = datetuple
    sec = round(sec)
    if sec == 60:
        return datetime.datetime(yr, mo, date, hr, mn, 59, tzinfo=tz_utc) + datetime.timedelta(seconds=1) 
    return datetime.datetime(yr, mo, date, hr, mn, sec, tzinfo=tz_utc)
    

def jday2datetime_us(jdt: float) -> datetime.datetime:
    """Turn julian day into datetime object to nearest microsecond"""
    datetuple = invjday(jdt)
    yr, mo, date, hr, mn, sec = datetuple
    sec, us = np.divmod(sec, 1)
    return datetime.datetime(yr, mo, date, hr, mn, int(sec), int(us*1e6), tzinfo=tz_utc)


def invjday_array(jd):
    """This procedure finds the year, month, day, hour, minute and second given the
    julian date. jd can be ut1, tdt, tdb, etc.
    Args:
        jd : float, julian date, days from 4713 BCE
    Outputs:
        year : int, year between 1900 - 2100
        mon : int, month between 1 - 12
        day : int, day between 1 - 31
        hr : int, hour between 0 - 23
        min : int, minute between 0 - 59
        sec : float, second between 0.0 - 59.999
    References:
        Vallado, 2007, 208, alg 22, ex 3-13
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # find year and days of the year
    jd = np.atleast_1d(jd)
    temp = jd - 2415019.5
    tu = temp / 365.25  # julian centuries from 0 h jan 0, 1900
    year = (1900 + np.floor_divide(tu, 1.0).astype(int)).astype(int)
    leapyrs = np.floor_divide(((year - 1901) * 0.25), 1.0).astype(int)  # number of leap years from 1900
    # optional nudge by 8.64x10-7 sec to get even outputs
    # day of year plus fractional portion of a day
    days = temp - ((year - 1900) * 365.0 + leapyrs) + 0.00000000001
    # check for case of beginning of a year
    day1_idx = days < 1.0
    year[day1_idx] = year[day1_idx] - 1
    leapyrs[day1_idx] = np.floor_divide((year[day1_idx] - 1901) * 0.25, 1.0).astype(int)
    days[day1_idx] = temp[day1_idx] - ((year[day1_idx] - 1900) * 365.0 + leapyrs[day1_idx])
    # find remaing data
    jd_size = len(jd)
    mon = np.empty(jd_size, dtype=int)
    day = np.empty(jd_size, dtype=int)
    hr = np.empty(jd_size, dtype=int)
    minute = np.empty(jd_size, dtype=int)
    sec = np.empty(jd_size)
    for i in range(jd_size):
        mon[i], day[i], hr[i], minute[i], sec[i] = days2mdhms(year[i], days[i])
    sec = sec - 0.00000086400
    return year, mon, day, hr, minute, sec


def jday2datetime_array(jdt_array):
    """Turn julian day into datetime object rounded to the second"""
    datetuple = invjday_array(jdt_array)
    yr, mo, date, hr, mn, sec = datetuple
    sec = np.around(sec).astype(int)
    dt_array = np.empty(yr.size, dtype=object)
    for i, data in enumerate(zip(yr, mo, date, hr, mn, sec)):
        dt_array[i] = datetime.datetime(*data, tzinfo=tz_utc)
    return dt_array


def jday2datetime_us_array(jdt_array):
    """Turn julian day into datetime object rounded to microsecond"""
    datetuple = invjday_array(jdt_array)
    yr, mo, date, hr, mn, sec = datetuple
    sec, us = np.divmod(sec, 1)
    sec = sec.astype(int)
    us *= 10**6
    dt_array = np.empty(yr.size, dtype=object)
    for i, data in enumerate(zip(yr, mo, date, hr, mn, sec, us.astype(int))):
        dt_array[i] = datetime.datetime(*data, tzinfo=tz_utc)
    return dt_array


def compute_time_array(dt_start: datetime.datetime, dt_end: datetime.datetime, dt_seconds: float) -> Time:
    """
    Create astropy Time object
    """
    jdt0 = julian_date(dt_start)
    jdtf = julian_date(dt_end)
    dt_days = dt_seconds/(24*60*60.0)
    jd_array = np.arange(jdt0, jdtf, dt_days, dtype=float)
    return Time(jd_array, format='jd')


def compute_time_array_from_date(date_start: datetime.date, date_end: datetime.date, dt_seconds: float) -> Time:
    """
    Create astropy Time object from python date
    """
    jd_array = julian_date_array_from_date(date_start, date_end, dt_seconds)
    return Time(jd_array, format='jd')


def julian_date_array_from_date(date_start: datetime.date, date_end: datetime.date, dt_seconds: float) -> np.ndarray:
    """
    Create a numpy array of julian date values from python dates
    """
    jdt0 = julian_day(date_start.year, date_start.month, date_start.day)
    jdtf = julian_day(date_end.year, date_end.month, date_end.day)
    dt_days = dt_seconds/(24*60*60.0)
    return np.arange(jdt0, jdtf, dt_days, dtype=float)