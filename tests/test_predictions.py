# Test using pytest

import datetime

from numpy.testing import assert_allclose, assert_almost_equal
import numpy as np
import pytest

from passpredict import predictions
from passpredict import timefn
from passpredict.constants import ASEC2RAD
from passpredict.utils import get_TLE
from passpredict.schemas import Location, Satellite
from passpredict.propagate import propagate_satellite




def test_satellite_visible():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([[-2811.2769, 3486.2632, 5069.5763]]).T  # ECI coords
    rsite = np.array([[-3414.0283, 3258.1636, 4276.1212]]).T  # ECI coords
    rho = np.array([[-773.8654, -581.4980, 328.8145]]).T  # SEZ coords
    dt = datetime.datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = np.array([timefn.julian_date2(dt)])
    vis = predictions.satellite_visible(rsat, rsite, rho, jdt)
    assert vis[0] > 2


@pytest.mark.predict
def test_predict():
    """
    Just confirm that the predict() function doesn't error
    """
    satellite = Satellite(id=25544, name='ISS')
    location = Location(lat=30.2711, lon=-97.7434, h=0, name='Austin, Texas')
    tle = get_TLE(satellite)
    dt_start = timefn.truncate_datetime(datetime.datetime.now())# - datetime.timedelta(days=1)
    dt_end = dt_start + datetime.timedelta(days=10)
    min_elevation = 10.01 # degrees
    overpasses = predictions.predict(location, satellite, dt_start=dt_start, dt_end=dt_end, dt_seconds=1, min_elevation=min_elevation, verbose=True)
    assert True


# Create a prediction test suite comparing results to 
#    1. heavens-above.com
#    2. calsky.com
#    3. n2Y0.com
#    4. skyfield
#    5. pyephem
#
# Use multiple satellites and multiple locations at different latitudes
#   Satellites: 
#     ISS (25544)
#     Hubble
#     Starlink-8
#     Lightsail 2
#     Envisat
#     Terra 
#     a geosynchronus sat
#     a retrograde orbit sat
#     a molniya orbit sat
#     a sun-synchronus orbit sat
#     a geostationary orbit sat
#
#   Locations (lat, lon):
#     Longyearbyen, Norway (78.2)
#     Inuvik, Canada (68.3607 , -133.7230)
#     Helsinki, Finland (60.17)
#     Kiev, Ukraine (50.45)
#     New York, NY, USA (40.7128, -74.0060)
#     Austin, Texas, USA (30.2672, -97.7431)
#     Mexico City, Mexico (~20)
#     Bissau, Guinea Bissau (~10)
#     Quito, Ecuador (~0)
#     Johannesburg, South Africa (-26)
#     Sydney, Australia (-33)
#     Christchurch, New Zealand (-43)
#     Punta Arenas, Chile (-53)
#     McMurdo Station, Antarctica (-77)
predict_results = {

}