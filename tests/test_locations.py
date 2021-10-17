# test locations.py
from numpy.testing import assert_allclose
import pytest

from passpredict.locations import Location


@pytest.mark.parametrize(
    'lat, lon, h, recef, tol',
    (
        # pytest.param(42.38, -71.38, 24, [1526.122, -4465.064, 4276.894], 1e-3, id='Vallado, Eg 11-6, p.912'),
        pytest.param(39.007, -104.883, 2187, [-1275.1219, -4797.9890, 3994.2975], 1e-4, id='Vallado, Eg 7-1, p.431'),
    )
)
def test_location_instance_ecef_position(lat, lon, h, recef, tol):
    location = Location("", latitude_deg=lat, longitude_deg=lon, elevation_m=h)
    assert_allclose(location.position_ecef, recef, atol=tol)
