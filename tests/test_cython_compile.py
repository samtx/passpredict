import pytest

from passpredict import Location, Satellite


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