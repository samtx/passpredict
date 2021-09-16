import pytest

from passpredict import Location, Satellite, OMM
from passpredict import predict
from passpredict import timefn


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
    ('tle1', 'tle2'),
    (
        (
            '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999',
            '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675',
        ),
    )
)
def test_satellite_cython_object_init_using_tle(tle1, tle2):
    """
    Create satellite object from orbit OMM data
    """
    omm = OMM.from_tle(tle1, tle2)
    satellite = Satellite(omm)
    assert satellite.inclo == omm.inclo
    assert satellite.ecco == omm.ecco
    assert satellite.nodeo == omm.nodeo
    assert satellite.argpo == omm.argpo
    assert satellite.mo == omm.mo
    assert satellite.bstar == omm.bstar
    assert satellite.jdsatepoch == omm.jdsatepoch
    assert satellite.jdsatepochF == omm.jdsatepochF
    assert satellite.no_kozai == omm.no_kozai
    assert satellite.elnum == omm.elnum


def test_compute_elevation_angle():
    """
    Test cython function works properly. Don't test for accuracy
    """
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    omm = OMM.from_tle(tle1, tle2)
    satellite = Satellite(omm)
    location = Location(32.1, -97.5, 20)
    jd = 2458871.5
    el = predict.compute_elevation_angle(jd, location, satellite)
    assert el >= -90
    assert el <= 90
    assert isinstance(el, float)


# @pytest.mark.parametrize(
#     ('tle1', 'tle2', 'brightness'),
#     (
#         (
#             '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999',
#             '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675',
#             -1.2,
#         ),
#     )
# )
# def test_satellite_cython_object_init_using_tle_set_brightness(tle1, tle2, brightness):
#     """
#     Create satellite object from orbit OMM data
#     """
#     omm = OMM.from_tle(tle1, tle2)
#     satellite = Satellite(omm)
#     satellite.brightness = brightness
#     assert satellite.brightness == brightness
