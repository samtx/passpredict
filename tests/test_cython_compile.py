import datetime

import pytest

from passpredict import Location, Satellite
import passpredict.timefn as timefn



@pytest.mark.parametrize(
    ('lat', 'lon', 'h', 'name'),
    (
        [32.1234, -97.124, 39, 'Austin'],
        [55.147, -45, 0, ''],
        [55.147, -45, 0, None],
    )

)
def test_location_cython_object_init(lat, lon, h, name):
    """
    Create a Location object
    """
    if name is not None:
        location = Location(lat, lon, h, name)
    else:
        location = Location(lat, lon, h)
    assert location.lat == lat
    assert location.lon == lon
    assert location.h == h
    assert location.name == name


@pytest.mark.parametrize(
    ('year', 'month', 'day', 'hr', 'minute', 'sec', 'jd_expected', 'jdfrac_expected'),
    (
        (2020, 1, 23, 3, 9, 45.4687, 2458871.5, 0.13177625810185184),
        (2020, 1, 23, 3, 9, 45, 2458871.5, 0.13177083333333334),
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



