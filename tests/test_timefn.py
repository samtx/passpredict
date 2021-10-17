# test _time.pyx
from datetime import datetime, timezone

import numpy as np
from numpy.testing import assert_allclose
import pytest

from passpredict._time import julian_date, epoch_to_jd, jday2datetime


def test_epoch_to_jd():
    """
    Vallado, p. 107
    """
    yr = 93
    dt = 352.53502934
    jd_actual, jdfr_actual = julian_date(1993, 12, 18, 12, 50, 26.5350)
    jd, jdfr = epoch_to_jd(yr, dt)
    assert_allclose(jd, jd_actual)
    assert_allclose(jdfr, jdfr_actual)


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
        assert_allclose(jd, jd_expected, atol=atol)

    def test_jday2datetime(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        """Convert a Julian date to a datetime and back"""
        dt_computed = jday2datetime(jd_expected)
        sec, us = np.divmod(sec, 1)
        dt_desired = datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        dt_difference = dt_computed - dt_desired
        assert abs(dt_difference.total_seconds()) < 0.5

    def test_julian_day_to_datetime_and_back(self, yr, mo, dy, hr, mn, sec, jd_expected, atol):
        jd, jdfr = julian_date(yr, mo, dy, hr, mn, sec)
        jd += jdfr
        dt = jday2datetime(jd)
        sec, us = np.divmod(sec, 1)
        delta = dt - datetime(yr, mo, dy, hr, mn, int(sec), int(us*1e6), tzinfo=timezone.utc)
        assert abs(delta.total_seconds()) < 0.5


if __name__ == "__main__":
    pytest.main(['-v', __file__])