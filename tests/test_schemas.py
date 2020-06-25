# Test passpredict/schemas.py

import passpredict.schemas as schemas
from passpredict.propagate import epoch_from_tle
import datetime
import pytest


def test_Satellite():
    satid = 25544
    name = "International Space Station"
    satellite = schemas.Satellite(id=satid, name=name)
    assert satellite.id == satid
    assert satellite.name == name


@pytest.mark.xfail(strict=True)
def test_Satellite_fail():
    """
    This test is expected to fail
    """
    satid = 25544
    name = "International Space Station"
    satellite = schemas.Satellite(satid, name)
    assert satellite.id == satid
    assert satellite.name == name


def test_Tle():
    """
    From Celestrak

    ISS (ZARYA)             
    1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999
    2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675
    """
    satid = 25544
    name = "International Space Station"
    satellite = schemas.Satellite(id=satid, name=name)
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    epoch = epoch_from_tle(tle1)
    tle = schemas.Tle(tle1=tle1, tle2=tle2, epoch=epoch, satellite=satellite)
    assert tle.tle1 == tle1
    assert tle.tle2 == tle2
    assert tle.epoch == epoch
    assert tle.satellite == satellite


def test_Point():
    dt = datetime.datetime(2020, 6, 1, 12, 34, 56, tzinfo=datetime.timezone.utc)
    azimuth = 270.5
    elevation = 84.1
    range_ = 3200.214
    point = schemas.Point(datetime=dt, azimuth=azimuth, elevation=elevation, range=range_)
    assert point.datetime == dt
    assert point.azimuth == azimuth
    assert point.elevation == elevation
    assert point.range == range_


def test_Overpass():
    start_pt = schemas.Point(
        datetime = datetime.datetime(2020, 6, 1, 12, 34, 56, tzinfo=datetime.timezone.utc),
        azimuth = 270.5,
        elevation = 12,
        range = 3200.214
    )
    max_pt = schemas.Point(
        datetime = datetime.datetime(2020, 6, 1, 12, 35, 59, tzinfo=datetime.timezone.utc),
        azimuth = 358,
        elevation = 84.1,
        range = 2500
    )
    end_pt = schemas.Point(
        datetime = datetime.datetime(2020, 6, 1, 12, 37, 12, tzinfo=datetime.timezone.utc),
        azimuth = 20,
        elevation = 10.0000001,
        range = 3987.1483321548
    )
    overpass = schemas.Overpass(
        start_pt=start_pt,
        max_pt=max_pt,
        end_pt=end_pt
    )
    assert overpass.start_pt == start_pt
    assert overpass.start_pt.datetime == datetime.datetime(2020, 6, 1, 12, 34, 56, tzinfo=datetime.timezone.utc)
    assert overpass.max_pt == max_pt


def test_Location():
    lat = 30.2672
    lon = -97.7431
    h = 88
    name = 'Austin, Texas'
    austin = schemas.Location(lat=lat, lon=lon, h=h, name=name)
    assert austin.lat == lat
    assert austin.lon == lon
    assert austin.h == h
    assert austin.name == name


if __name__ == "__main__":
    pytest.main(['-v', __file__])
