# test timefn.py
from passpredict import timefn
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal
from datetime import datetime, timedelta
import pytz


def test_utc2tt():
    """Vallado, Eg. 3-7"""
    # Mountain standard time (UTC-6)
    dt = datetime(2004, 5, 14, 10, 43, 0)
    deltaAT = 32  # sec
    deltaUT1 = -0.463326  # sec
    # get UTC by adjusting from MST timezone
    dt += timedelta(hours=6)
    jd_utc = timefn.julian_date(dt)
    tt = timefn.utc2tt(jd_utc, deltaAT=deltaAT, deltaUT1=deltaUT1)
    assert_almost_equal(tt, 0.043674121031, decimal=12)


def test_julian_date():
    """
    Vallado, eg.3-4
    """
    yr, mo, dy = 1996, 10.0, 26.0
    hr, mn, sec = 14.0, 20.0, 0.0
    jd = timefn.julian_date(yr, mo, dy, hr, mn, sec)
    jdT = 2450383.09722222
    assert_almost_equal(jd, jdT, decimal=8)


def test_julian_date_datetime():
    """
    Vallado, eg.3-4
    """
    yr, mo, dy = 1996, 10, 26
    hr, mn, sec = 14, 20, 0
    dt = datetime(yr, mo, dy, hr, mn, sec)
    jd = timefn.julian_date(dt)
    jdT = 2450383.09722222
    assert_almost_equal(jd, jdT, decimal=8)


# def test_julian_date_datetime2():
#     """
#     Vallado, eg. 3-15, p. 230
#     """
#     dt = datetime(2004, 4, 6, 7, 51, )
#     jd = timefn.julian_date(dt)
#     jdT = 2450383.09722222
#     assert_almost_equal(jd, jdT, decimal=8)


def test_julian_date_vectorized():
    """Use an array of datetimes to find the Julian Date"""
    dt_ary = np.arange(
        "2019-09-14T00:00:00", "2019-10-07T00:00:00", 200, dtype="datetime64"
    )
    jd_vectorized = np.vectorize(timefn.julian_date)
    jd_ary = jd_vectorized(dt_ary)
    # print(jd_ary)


def test_jd_from_skyfield():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)
    assert_almost_equal(jd, 2453101.8274067827, decimal=12)


def test_jd_from_skyfield2():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date2(dt)
    jd_desired = 2453101.8274067827
    print(f"jd   ={jd:20.15f}")
    print(f"jddes={jd_desired:20.15f}")
    print(f"diff ={jd-jd_desired:0.15f}")
    assert_almost_equal(jd, 2453101.8274067827, decimal=12)


def test_jd_from_skyfield3():
    """From skyfield.tests.test_earth_satellites.py"""
    sec = 28.386 - 0.439961
    yr, mo, dy = 2004, 4, 6
    hr, mn = 7, 51
    jd = timefn.julian_date2(yr, mo, dy, hr, mn, sec)
    jd_desired = 2453101.8274067827
    print(f"jd   ={jd:20.15f}")
    print(f"jddes={jd_desired:20.15f}")
    print(f"diff ={jd-jd_desired:0.15f}")
    assert_almost_equal(jd, jd_desired, decimal=12)

def test_jday2datetime():
    """Convert a Julian date to a datetime and back"""
    from datetime import timezone
    # jd = timefn.julian_date(dt)
    julian_date = 2450383.09722222  # 1996-10-26 14:20:0
    dt_computed = timefn.jday2datetime(julian_date)
    dt_desired = datetime(1996, 10, 26, 14, 20, 0, tzinfo=timezone.utc)
    dt_difference = dt_computed - dt_desired
    assert abs(dt_difference.total_seconds()) < 500e-6  # 500 microseconds

if __name__ == "__main__":
    # test_jd_from_skyfield()
    # test_jd_from_skyfield2()
    # test_jd_from_skyfield3()
    test_jday2datetime()