import pytest
from pytest import approx

from passpredict import geocoding


def test_geocoder_nominatium():
    query_str = "austin, texas"
    geocoder = geocoding.NominatimGeocoder()
    location = geocoder.query(query_str)
    assert location.latitude_deg == approx(30.2711)
    assert location.longitude_deg == approx(-97.7437)
    assert "Austin" in location.name
    assert "Texas" in location.name


if __name__ == "__main__":
    pytest.main([__file__, "-v"])