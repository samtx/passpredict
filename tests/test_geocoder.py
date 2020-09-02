import pytest

from passpredict import geocoding


def test_geocoder_nominatium():
    query_str = "austin, texas"
    res = geocoding.geocoder(query_str)
    lat = round(float(res['lat']), 4)
    lon = round(float(res['lon']), 4)
    assert 30.2 <= lat <= 30.3
    assert -97.7 >= lon >= -97.8