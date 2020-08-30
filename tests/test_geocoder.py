import pytest

from passpredict import geocoding


def test_geocoder_nominatium():
    query_str = "austin, texas"
    res = geocoding.geocoder(query_str)
    lat = round(float(res['lat']), 4)
    lon = round(float(res['lon']), 4)
    assert lat == 30.2672
    assert lon == -97.7431