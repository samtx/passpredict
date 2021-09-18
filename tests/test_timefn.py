import datetime

import pytest

from passpredict import timefn


@pytest.mark.parametrize(
    ('year', 'month', 'day', 'hr', 'minute', 'sec', 'jd_expected', 'jdfrac_expected'),
    (
        (2020, 1, 23, 3, 9, 45.4687, 2458871.5, 0.13177625810185184),
        (2020, 1, 23, 3, 9, 45, 2458871.5, 0.13177083333333334),
        (1997, 4,  2, 1, 8,  0, 2450540.5, 0.047222222),
    )
)
def test_julian_date(year, month, day, hr, minute, sec, jd_expected, jdfrac_expected):
    """
    Compute the julian date
    """
    jd, jdfrac = timefn.julian_date(year, month, day, hr, minute, sec)
    print(jd, jdfrac)
    assert jd == pytest.approx(jd_expected)
    assert jdfrac == pytest.approx(jdfrac_expected)


@pytest.mark.parametrize(
    ('jd', 'datetime_expected'),
    (
        (2458871.5 + 0.13177625810185184, datetime.datetime(2020, 1, 23, 3, 9, 45, tzinfo=datetime.timezone.utc)),
    )
)
def test_jday2datetime_cython(jd, datetime_expected):
    """
    Create a datetime object from julian date float
    Datetime object is rounded to nearest second
    """
    datetime_ = timefn.jday2datetime(jd)
    assert datetime_ == datetime_expected