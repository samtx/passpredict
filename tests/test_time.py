# test time.py and _time.pyx
from datetime import datetime, timezone

import numpy as np
import pytest
from pytest import approx

from passpredict import time as ptime
from passpredict import _time as _ptime
from passpredict.tle import jd_to_epoch_string
from passpredict.constants import DAYSEC

def test_epoch_to_jd():
    """
    Vallado, p. 107
    """
    yr = 93
    dt = 352.53502934
    jd_actual, jdfr_actual = _ptime.julian_date(1993, 12, 18, 12, 50, 26.53502934)
    jd, jdfr = _ptime.epoch_to_jd(yr, dt)
    assert jd == approx(jd_actual)
    assert jdfr == approx(jdfr_actual)


def test_jd_to_epoch_string():
    """
    Vallado, p. 107
    """
    jd = 2449340.03502934
    epoch_str = jd_to_epoch_string(jd)
    assert epoch_str == "93352.53502934"


@pytest.mark.parametrize(
    'yr, mo, dy, hr, mn, sec, jd_expected, atol',
    [
        pytest.param(1996, 10, 26, 14, 20,         0, 2450383.09722222, 1e-8, id='Vallado, eg.3-4'),
        pytest.param(2004,  4,  6,  7, 52, 32.570009, 2453101.828154745, 1e-9, id='Vallado, eg.3-15'),
        pytest.param(2004,  4,  6,  7, 51, 27.946039, 2453101.8274067827, 1e-10, id='skyfield.tests.test_earth_satellites'),
    ]
)
class TestTimeFunctions:

    def test_julian_date_from_components(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        jd, jdfr = _ptime.julian_date(yr, mo, dy, hr, mn, sec)
        jd += jdfr
        assert jd == approx(jd_expected, abs=atol)

    def test_jday2datetime(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        """Convert a Julian date to a datetime and back"""
        dt_computed = _ptime.jday2datetime(jd_expected)
        sec, us = divmod(sec, 1)
        dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        dt_difference = dt_computed - dt_desired
        assert dt_difference.total_seconds() == approx(0.0, abs=0.5)

    def test_jday2datetime_us(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        """Convert a Julian date to a datetime with microseconds and back"""
        dt_computed = _ptime.jday2datetime_us(jd_expected)
        sec, us = divmod(sec, 1)
        dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        dt_difference = dt_computed - dt_desired
        assert dt_difference.total_seconds() == approx(0.0, abs=0.0005)

    def test_julian_day_to_datetime_and_back(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        jd, jdfr = _ptime.julian_date(yr, mo, dy, hr, mn, sec)
        jd += jdfr
        dt = _ptime.jday2datetime(jd)
        sec, us = divmod(sec, 1)
        delta = dt - datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        assert delta.total_seconds() == approx(0.0, abs=0.5)

    def test_julian_day_to_datetime_us_and_back(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        jd, jdfr = _ptime.julian_date(yr, mo, dy, hr, mn, sec)
        jd += jdfr
        dt = _ptime.jday2datetime_us(jd)
        sec, us = divmod(sec, 1)
        delta = dt - datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        assert delta.total_seconds() == approx(0.0, abs=0.5)

@pytest.mark.parametrize(
    'jd, jd_expected',
    (
        (2450383.09722222, 2450383.09722222),
        (2453101.828154745, 2453101.828148148),
        (2453101.8274067827, 2453101.8273958336),
    )
)
def test_julian_date_round_to_second(jd, jd_expected):
    jd2 = ptime.julian_date_round_to_second(jd)
    tol = 1/DAYSEC #  1/86400
    assert jd2 == approx(jd_expected, abs=tol)
    _, rem = divmod(jd2, tol)
    assert rem < tol



@pytest.mark.parametrize('mjd, jd, jdfr', [
    (59616.42243220018, 2459617.0, -0.07756779981481476),
    (59616.50286387156, 2459617.0, 0.002863871562500009),
    (50382.59722222015, 2450383.09722222, 0),
    (53101.328154745046, 2453101.828154745, 0),
    (53101.32740678266, 2453101.8274067827, 0),
    (59617.054584761194, 2459617.554584761, 0)
])
def test_mjd2jdfr(mjd, jd, jdfr):
    _jd, _jdfr = _ptime.mjd2jdfr(mjd)
    _jdsum = _jd + _jdfr
    jdsum = jd + jdfr
    assert _jdsum == approx(jdsum, abs=1e-12)
    assert abs(_jdfr) < 1


@pytest.mark.parametrize(
    'yr, mo, dy, hr, mn, sec, jd_expected, mjd_expected, atol',
    [
        pytest.param(1996, 10, 26, 14, 20,         0, 2450383.09722222, 50382.59722222015, 1e-8, id='Vallado, eg.3-4'),
        pytest.param(2004,  4,  6,  7, 52, 32.570009, 2453101.828154745, 53101.328154745046, 1e-9, id='Vallado, eg.3-15'),
        pytest.param(2004,  4,  6,  7, 51, 27.946039, 2453101.8274067827, 53101.32740678266, 1e-10, id='skyfield.tests.test_earth_satellites'),
    ]
)
class TestMjd2Datetime:
    def test_mjd2datetime(self, yr, mo, dy, hr, mn, sec, jd_expected, mjd_expected, atol):
        """Convert a modified julian date to a datetime"""
        dt_computed = _ptime.mjd2datetime(mjd_expected)
        sec, us = divmod(sec, 1)
        dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        dt_difference = dt_computed - dt_desired
        assert dt_difference.total_seconds() == approx(0.0, abs=0.0005)

    def test_mjd2datetime_and_back(self, yr, mo, dy, hr, mn, sec, jd_expected, mjd_expected, atol):
        """Convert a modified julian date to a datetime with microseconds and back"""
        dt1 = _ptime.mjd2datetime(mjd_expected)
        mjd = _ptime.datetime2mjd(dt1)
        dt2 = _ptime.mjd2datetime(mjd)
        dt_difference = dt1 - dt2
        assert dt_difference.total_seconds() == approx(0.0, abs=1e-6)

    def test_datetime2mjd(self, yr, mo, dy, hr, mn, sec, jd_expected, mjd_expected, atol):
        """  Convert a datetime to modified julian date  """
        sec, us = divmod(sec, 1)
        dt = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        mjd = _ptime.datetime2mjd(dt)
        assert atol < 9.999e-7  # ensure tolerance is at least 1e-8
        assert mjd == approx(mjd_expected, abs=atol)


if __name__ == "__main__":
    pytest.main(['-v', __file__])