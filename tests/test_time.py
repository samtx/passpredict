# test time.py and _time.pyx
from datetime import datetime, timezone
from passpredict.constants import DAYSEC

import numpy as np
import pytest
from pytest import approx

from passpredict import *
import passpredict.time
from passpredict._time import jday2datetime_us
from passpredict.tle import jd_to_epoch_string

def test_epoch_to_jd():
    """
    Vallado, p. 107
    """
    yr = 93
    dt = 352.53502934
    jd_actual, jdfr_actual = julian_date(1993, 12, 18, 12, 50, 26.53502934)
    jd, jdfr = epoch_to_jd(yr, dt)
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
        pytest.param(1996, 10, 26, 14, 20,         0, 2450383.09722222,   1e-8, id='Vallado, eg.3-4'),
        pytest.param(2004,  4,  6,  7, 52, 32.570009, 2453101.828154745,  1e-9, id='Vallado, eg.3-15'),
        pytest.param(2004,  4,  6,  7, 51, 27.946039, 2453101.8274067827, 1e-10, id='skyfield.tests.test_earth_satellites'),
    ]
)
class TestTimeFunctions:

    def test_julian_date_from_components(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        jd, jdfr = julian_date(yr, mo, dy, hr, mn, sec)
        jd += jdfr
        assert jd == approx(jd_expected, abs=atol)

    def test_jday2datetime(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        """Convert a Julian date to a datetime and back"""
        dt_computed = jday2datetime(jd_expected)
        sec, us = np.divmod(sec, 1)
        dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        dt_difference = dt_computed - dt_desired
        assert dt_difference.total_seconds() == approx(0.0, abs=0.5)

    def test_jday2datetime_us(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        """Convert a Julian date to a datetime with microseconds and back"""
        dt_computed = jday2datetime_us(jd_expected)
        sec, us = np.divmod(sec, 1)
        dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        dt_difference = dt_computed - dt_desired
        assert dt_difference.total_seconds() == approx(0.0, abs=0.005)

    def test_julian_day_to_datetime_and_back(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        jd, jdfr = julian_date(yr, mo, dy, hr, mn, sec)
        jd += jdfr
        dt = jday2datetime(jd)
        sec, us = np.divmod(sec, 1)
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
    jd2 = passpredict.time.julian_date_round_to_second(jd)
    tol = 1/DAYSEC #  1/86400
    assert jd2 == approx(jd_expected, abs=tol)
    _, rem = divmod(jd2, tol)
    assert rem < tol


if __name__ == "__main__":
    pytest.main(['-v', __file__])