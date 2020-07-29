import pickle
from datetime import datetime, timedelta
from typing import List

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
from astropy import units as u
from astropy.coordinates import TEME, CartesianRepresentation, ITRS, get_sun
from astropy.time import Time

from .rotations import ecef2sez
from .solar import sun_pos, is_sat_illuminated, compute_sun_data
from .propagate import propagate_satellite, compute_satellite_data
from .timefn import julian_date, compute_time_array
from .schemas import Overpass, Satellite, Location, Tle
from .models import Sun, RhoVector, Sat
from .utils import get_TLE


def find_overpasses(location: Location, sats: List[Sat], times: Time, sun: Sun, min_elevation: float = 10) -> List[Overpass]:
    """
    Real-time computation for finding satellite overpasses of a topographic location.
    Can support multiple satellites over a single location
    """
    store_sat_id = True if len(sats) > 0 else False
    overpasses = []
    for sat in sats:
        rho = RhoVector(sat, location, sun)
        sat_overpasses = rho.find_overpasses(min_elevation, store_sat_id)
        overpasses += sat_overpasses
    return overpasses


def predict(location, satellite, dt_start=None, dt_end=None, dt_seconds=1, min_elevation=None, tle=None, cache=None, verbose=False, store_sat_id=False):
    """
    Full prediction algorithm:
      1. Download TLE data
      2. Propagate satellite using SGP4
      3. Predict overpasses based on site location
      4. Return overpass object and print to screen

    Params:
        location : Location object
            latitude of site location, in decimal, north is positive
        satellite: Satellite object
            satellite ID number in Celestrak, ISS is 25544
    """
    if dt_start is None:
        dt_start = datetime.now()
    if dt_end is None:
        dt_end = dt_start + timedelta(days=14)
    if tle is None:
        tle = get_TLE(satellite)

    t = compute_time_array(dt_start, dt_end, dt_seconds)
    
    if verbose:
        print("Compute sun position...")
    sun = compute_sun_data(t)

    if verbose:
        print(f"Begin propagation from {dt_start.isoformat()} to {dt_end.isoformat()}...")
    sat = compute_satellite_data(tle, t, sun)

    if verbose:
        print("Compute sun-satellite illumination...")
    sat.illuminated = is_sat_illuminated(sat.rECEF, sun.rECEF)

    if verbose:
        print('begin prediction...')
    overpasses = find_overpasses(location, [sat], t, sun, min_elevation)
    return overpasses
    


