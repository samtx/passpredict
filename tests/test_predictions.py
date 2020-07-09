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
    satellite = Satellite(id=25544, name='ISS')
    location = Location(lat=30.2711, lon=-97.7434, h=0, name='Austin, Texas')
    tle = get_TLE(satellite)
    dt_start = timefn.truncate_datetime(datetime.datetime.now())# - datetime.timedelta(days=1)
    dt_end = dt_start + datetime.timedelta(days=10)
    min_elevation = 10.01 # degrees
    overpasses = predictions.predict(location, satellite, dt_start=dt_start, dt_end=dt_end, dt_seconds=1, min_elevation=min_elevation, verbose=True)
    assert True


# def test_predict_passes():
#     """
#     Right now, just test that it doesn't error
#     """
#     # Set up satellite position
#     dt_seconds = 10
#     num_days = 5
#     min_elevation = 10.0

#     austin = Location(lat=30.2672, lon=-97.7431, h=0, name='Austin')
#     iss = Satellite(id=25544, name='ISS')
#     iss_tle = get_TLE(iss)
#     datetime_start = timefn.truncate_datetime(
#         datetime.datetime.now(tz=datetime.timezone.utc)
#     )
#     datetime_end = datetime_start + datetime.timedelta(days=num_days)
#     iss_rv = propagate_satellite.__wrapped__(
#         iss_tle.tle1, iss_tle.tle2, datetime_start, datetime_end, dt_seconds
#     )
#     overpasses = predictions.predict_passes(
#         austin.lat, austin.lon, austin.h,
#         iss_rv.rECEF, iss_rv.rECI, iss_rv.julian_date,
#         min_elevation=min_elevation
#     )
#     assert True



# def test_riseset():
#     """Eg. ,
#     Test orbit propagation
#     """
#         # Obj, n [rev/solar day],         e,  i [deg],  w, Omega, M
#     orbital_objects = [
#         (   1,        1.00272141, 0.0000032,  00.0956, 0.,    0., 0.),
#         (   2,        8.36589235, 0.0080158,  90.0175, 0.,    0., 0.),
#         (   3,        0.24891961, 0.9363060,  64.9874, 0.,    0., 0.),
#         (   4,        0.21467209, 0.0668128,  57.3500, 0.,    0., 0.),
#         (   5,       13.37659679, 0.0145072,  90.2619, 0.,    0., 0.),
#         (   6,       16.09769232, 0.0078742,  82.8709, 0.,    0., 0.),
#         (   7,        1.00271920, 0.0003109,  00.0099, 0.,    0., 0.),
#         (   8,       12.41552416, 0.0036498,  74.0186, 0.,    0., 0.),
#         (   9,       13.84150848, 0.0048964, 144.6414, 0.,    0., 0.)
#     ]



