# Test using pytest

from datetime import datetime, timedelta, timezone

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
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
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
    dt_start = timefn.truncate_datetime(datetime.now())# - datetime.timedelta(days=1)
    dt_end = dt_start + timedelta(days=2)
    min_elevation = 10.01 # degrees
    overpasses = predictions.predict(location, satellite, dt_start=dt_start, dt_end=dt_end, dt_seconds=1, min_elevation=min_elevation, verbose=True)
    assert True


# @pytest.mark.predict
# def test_predict():
#     """
#     Use test case in Vallado Fortran code 
        
#     TESTASTF.FOR, line 3066, TESTPREDICT
#     ----
#     SUBROUTINE TESTPREDICT
#         IMPLICIT NONE
#         REAL*8 JD, Latgd, LST, r(3),v(3),RS(3), Rho, Az, El, trtasc,
#      &         tdecl
#         CHARACTER WhichKind
#         CHARACTER*11 Vis
#         INCLUDE 'astmath.cmn'

#         Read(10,*) JD, Latgd, LST,  WhichKind
#         Write(20,*) JD, Latgd, LST, WhichKind
#         Read(10,*) r(1), r(2), r(3)
#         Write(20,*) r(1), r(2), r(3)
#         Read(10,*) v(1), v(2), v(3)
#         Write(20,*) v(1), v(2), v(3)
#         Read(10,*) rs(1), rs(2), rs(3)
#         Write(20,*) rs(1), rs(2), rs(3)
#         Latgd = Latgd * deg2rad
#         LST = LST * deg2rad

#         CALL PREDICT       ( JD,latgd,LST, r,v,rs, WhichKind,
#      &                           Rho,Az,El,tRtasc,tDecl, Vis )

#         Write(20,*) '  Results:'
#         Write(20,*) Rho*6378.1363D0,Az*Rad2deg,El*rad2deg,tRtasc,
#      &              tDecl, Vis

#       RETURN
#       END ! Predict


#     TESTASTF.DAT, line 775, case 503
#     ----
#     503  xx Predict
#          2451524.2355765  35.23234  68.3294875 S
#             1.02     0.903    1.04
#             0.784    0.27    0.12
#             0.238756  0.592834 0.9092385
#     """
#     jd = 2451524.2355765
#     lat = 35.23234
#     lst = 68.3294875
#     r = np.array(((1.02), (0.903), (1.04)))  # position vector in km
#     v = np.array(((0.784), (0.27), (0.12)))  # velocity vector in km
#     rsite = np.array(((0.238756), (0.592834), (0.9092385)))  # site position vector in km

#     satellite = Satellite(id=25544, name='ISS')
#     location = Location(lat=30.2711, lon=-97.7434, h=0, name='Austin, Texas')
#     tle = get_TLE(satellite)
#     dt_start = timefn.truncate_datetime(datetime.datetime.now())# - datetime.timedelta(days=1)
#     dt_end = dt_start + datetime.timedelta(days=10)
#     min_elevation = 10.01 # degrees
#     overpasses = predictions.predict(location, satellite, dt_start=dt_start, dt_end=dt_end, dt_seconds=1, min_elevation=min_elevation, verbose=True)
#     assert True


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