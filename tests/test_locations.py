# test rotations.py

import numpy as np
import pytest
from pytest import approx

from passpredict import Location

try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo


@pytest.mark.parametrize(
    'lat, lon, h, ecef_expected, tol',
    (
        pytest.param(42.38, -71.13, 24, [1526.122, -4465.064, 4276.894], 1e-3, id="Vallado, Eg. 11-6, p.912"),
        pytest.param(39.007, -104.883, 2187, [-1275.1219, -4797.9890, 3994.2975], 1e-4, id="Vallado, Eg. 7-1, p.431"),
    )
)
def test_location_instance_ecef_position(lat, lon, h, ecef_expected, tol):
    """
    Test Location object to create ecef position vector
    """
    location = Location("", latitude_deg=lat, longitude_deg=lon, elevation_m=h)
    ecef = location.recef
    assert ecef == approx(ecef_expected, abs=tol)


@pytest.mark.parametrize(
    'name, lat, lon, zoneinfo_str, offset_tuple',
    (
        pytest.param('Austin', 30.2672, -97.7431, 'America/Chicago', (-5, -6,), id='Austin, Texas'),
        pytest.param('London', 51.5072, 0.1276, 'Europe/London', (0, +1,), id='London, UK'),
        pytest.param('Phoenix', 33.4484, -112.0740, 'America/Phoenix', (-7,), id='Phoenix, Arizona'),
        pytest.param('Petropavlovsk-Kamchatsky', 53.01667, 158.65, 'Asia/Kamchatka', (+12,), id="Petropavlovsk-Kamchatsky, Kamchatka"),
    )
)
def test_location_timezone(name, lat, lon, zoneinfo_str, offset_tuple):
    location = Location(name, lat, lon, 0)
    assert location.timezone == ZoneInfo(zoneinfo_str)
    assert location.offset in offset_tuple


