from datetime import datetime, timedelta
from typing import List, Callable
import time

from astropy.time import Time
from numpy import ndarray

from .solar import compute_sun_data
from .propagate import compute_satellite_data
from .timefn import julian_date_array_from_date
from .schemas import Overpass, Location, Satellite
from .models import Sun, RhoVector, Sat
from .tle import get_TLE


def find_overpasses(
    jd: ndarray, 
    location: Location, 
    sats: List[Sat], 
    sun: Sun, 
    min_elevation: float = 10.0, 
    visible_only: bool =False
) -> List[Overpass]:
    """
    Real-time computation for finding satellite overpasses of a topographic location.
    Can support multiple satellites over a single location
    """
    store_sat_id = True if len(sats) > 1 else False
    overpasses = []
    for sat in sats:
        rho = RhoVector(jd, sat, location, sun)
        sat_overpasses = rho.find_overpasses(min_elevation, store_sat_id, visible_only=visible_only)
        overpasses += sat_overpasses
    return overpasses


def predict(
    location: Location,
    satellite: Satellite,
    date_start: datetime,
    date_end: datetime,
    min_elevation: float = 10.0,
):
    """
    Full prediction algorithm:
        - Propagate satellite using SGP4
        - Predict overpasses based on site location
    
    Params:
        location : Location object
            latitude of site location, in decimal, north is positive
        satellite: Satellite object
            satellite ID number in Celestrak, ISS is 25544
    """
    t = Time(jd, format='jd')
    sun = compute_sun_data(t)
    tle = get_TLE(satellite.id)
    sat = compute_satellite_data(tle, t, sun)
    return overpasses



