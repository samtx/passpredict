# Cython file for time functions
import numpy as np
cimport numpy as np
import datetime

cdef  (int, int, int, int, double) _days2mdhms(int year, double days):
    """This procedure converts the day of the year, days, to the equivalent
    month, day, hour, minute and second.
    Args:
        year : int, year between 1900 - 2100
        days : float, julian day of the year between 0.0 - 366.0
    Outputs:
        mon : int, month, 1 .. 12
        day : int, day, 1 .. 28,29,30,31
        hr : int, hour, 0 .. 23
        minute : int, minute 0 .. 59
        sec : float, second, 0.0 .. 59.999
    References:
        Vallado
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # non vectorized version...
    cdef int[:] lmonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    cdef int i, inttemp
    cdef double temp
    cdef int mon, day, hr, minute
    cdef double sec
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
    return (mon, day, hr, minute, sec)

# cdef double _invjday(double julian_date):
#     # find year and days of the year
#     julian_date = np.atleast_1d(julian_date)
#     temp = julian_date - 2415019.5
#     tu = temp / 365.25  # julian centuries from 0 h jan 0, 1900
#     year = 1900 + np.floor_divide(tu, 1.0).astype(int)
#     leapyrs = np.floor_divide(((year - 1901) * 0.25), 1.0).astype(int)  # number of leap years from 1900
#     # optional nudge by 8.64x10-7 sec to get even outputs
#     # day of year plus fractional portion of a day
#     days = temp - ((year - 1900) * 365.0 + leapyrs) + 0.00000000001
#     # check for case of beginning of a year
#     day1_idx = days < 1.0
#     year[day1_idx] = year[day1_idx] - 1
#     leapyrs[day1_idx] = np.floor_divide((year[day1_idx] - 1901) * 0.25, 1.0).astype(int)
#     days[day1_idx] = temp[day1_idx] - ((year[day1_idx] - 1900) * 365.0 + leapyrs[day1_idx])
#     # find remaing data
#     jd_size = len(julian_date)
#     mon = np.empty(jd_size)
#     day = np.empty(jd_size)
#     hr = np.empty(jd_size)
#     minute = np.empty(jd_size)
#     sec = np.empty(jd_size)

def jday2datetime(double[:] julian_date):
    """This procedure finds the year, month, day, hour, minute and second given the
    julian date. julian_date can be ut1, tdt, tdb, etc.
    Args:
        julian_date : float (n)
            julian date, days from 4713 BCE
    Outputs:
        datetime_array: datetime.datetime (n)
            array of datetime.datetime objects
    References:
        Vallado, 2007, 208, alg 22, ex 3-13
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # find year and days of the year
    cdef int i, jd_size
    cdef int mon_i, day_i, hr_i, minute_i
    cdef double sec_i
    julian_date = np.atleast_1d(julian_date)
    temp = julian_date - 2415019.5
    tu = temp / 365.25  # julian centuries from 0 h jan 0, 1900
    year = 1900 + np.floor_divide(tu, 1.0).astype(int)
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
    jd_size = len(julian_date)
    mon = np.empty(jd_size)
    day = np.empty(jd_size)
    hr = np.empty(jd_size)
    minute = np.empty(jd_size)
    sec = np.empty(jd_size)
    for i in range(jd_size):
         (mon_i, day_i, hr_i, minute_i, sec_i) = _days2mdhms(year[i], days[i])
         mon[i] = mon_i
         day_i[i] = day_i
         hr_i[i] = hr_i
         minute[i] = minute_i
         sec[i] = sec_i
    sec = sec - 0.00000086400
    return year, mon, day, hr, minute, sec