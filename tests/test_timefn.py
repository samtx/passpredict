# test timefn.py
import numpy as np
import pytest

from passpredict import timefn
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal
from datetime import datetime, timedelta, timezone


jd_params = [
    pytest.param(1996, 10, 26, 14, 20,         0, 2450383.09722222,   8, id='Vallado, eg.3-4'),
    pytest.param(2004,  4,  6,  7, 52, 32.570009, 2453101.828154745,  9, id='Vallado, eg.3-15'),
    pytest.param(2004,  4,  6,  7, 51, 27.946039, 2453101.8274067827, 10, id='skyfield.tests.test_earth_satellites'),    
]


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


@pytest.mark.parametrize(
    'yr, mo, dy, hr, mn, sec, jd_expected, decimal', jd_params
)
def test_julian_date_from_components(yr, mo, dy, hr, mn, sec, jd_expected, decimal):
    jd = timefn.julian_date(yr, mo, dy, hr, mn, sec)
    assert_almost_equal(jd, jd_expected, decimal=decimal)


@pytest.mark.parametrize(
    'yr, mo, dy, hr, mn, sec, jd_expected, decimal', jd_params
)
def test_julian_date_datetime(yr, mo, dy, hr, mn, sec, jd_expected, decimal):
    sec, us = np.divmod(sec, 1.)
    dt = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6))
    jd = timefn.julian_date(dt)
    assert_almost_equal(jd, jd_expected, decimal=decimal)


def test_julian_date_vectorized():
    """Use an array of datetimes to find the Julian Date"""
    dt_ary = np.arange(
        "2019-09-14T00:00:00", "2019-10-07T00:00:00", 200, dtype="datetime64"
    )
    jd_vectorized = np.vectorize(timefn.julian_date)
    jd_ary = jd_vectorized(dt_ary)


@pytest.mark.parametrize(
    'yr, mo, dy, hr, mn, sec, jd, decimal', jd_params
)
def test_jday2datetime(yr, mo, dy, hr, mn, sec, jd, decimal):
    """Convert a Julian date to a datetime and back"""
    from datetime import timezone
    dt_computed = timefn.jday2datetime(jd)
    sec, us = np.divmod(sec, 1)
    dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
    dt_difference = dt_computed - dt_desired
    assert abs(dt_difference.total_seconds()) < 500e-6  # 500 microseconds


if __name__ == "__main__":
    import pytest
    pytest.main(['-v', __file__])