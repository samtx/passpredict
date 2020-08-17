import pickle
from datetime import datetime, timedelta
from typing import List
import time

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
from astropy import units as u
from astropy.coordinates import TEME, CartesianRepresentation, ITRS, get_sun
from astropy.time import Time

from .rotations import ecef2sez
from .solar import sun_pos, is_sat_illuminated, compute_sun_data
from .propagate import propagate_satellite, compute_satellite_data
from .timefn import julian_date, compute_time_array, julian_day, compute_time_array_from_date
from .schemas import Overpass, Satellite, Location, Tle
from .models import Sun, RhoVector, Sat, SatPredictData, SunPredictData
from .utils import get_TLE, Cache


def find_overpasses(t: Time, location: Location, sats: List[Sat], sun: Sun, min_elevation: float = 10) -> List[Overpass]:
    """
    Real-time computation for finding satellite overpasses of a topographic location.
    Can support multiple satellites over a single location
    """
    store_sat_id = True if len(sats) > 1 else False
    overpasses = []
    for sat in sats:
        rho = RhoVector(t, sat, location, sun)
        sat_overpasses = rho.find_overpasses(min_elevation, store_sat_id)
        overpasses += sat_overpasses
    return overpasses


def predict(location, satellite, date_start=None, date_end=None, dt_seconds=1, min_elevation=None, cache=None, verbose=False, store_sat_id=False, print_fn=print):
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
    if date_start is None:
        date_start = datetime.date.today()
    if date_end is None:
        date_end = date_start + timedelta(days=14)
    
    if cache is None:
        t = compute_time_array_from_date(date_start, date_end, dt_seconds)
        sun = compute_sun_data(t)
        tle = get_TLE(satellite.id)
        sat = compute_satellite_data(tle, t, sun)

    else:

        with cache:

            time_key = date_start.strftime('%Y%m%d') + date_end.strftime('%Y%m%d') + str(dt_seconds)
            
            t = cache.get(time_key)
            if t is None:
                t = compute_time_array_from_date(date_start, date_end, dt_seconds)
                cache.set(time_key, t, ttl=86400)
            
            sun_key = 'sun_' + time_key
            sun = cache.get(sun_key)
            if sun is None:
                if verbose:
                    print_fn("Compute sun position... ", end=' ')
                t0 = time.perf_counter()
                sun = compute_sun_data(t)
                tf = time.perf_counter() - t0
                if verbose:
                    print_fn(f'{tf:0.3f} sec')
                cache.set(sun_key, sun, ttl=86400)

            tle_key = str(satellite.id) + '_tle'
            tle = cache.get(tle_key)
            if tle is None:
                tle = get_TLE(satellite.id)
                cache.set(tle_key, tle, ttl=86400)

            sat_key = tle_key + time_key + '_sat'
            sat = cache.get(sat_key)
            if sat is None:    
                if verbose:
                    print_fn(f"Propagate satellite from {date_start.isoformat()} to {date_end.isoformat()}... ", end=' ')
                t0 = time.perf_counter()
                sat = compute_satellite_data(tle, t, sun)
                tf = time.perf_counter() - t0
                if verbose:
                    print_fn(f'{tf:0.3f} sec')
                cache.set(sat_key, sat, ttl=86400)

    if verbose:
        print_fn('Predict overpasses... ', end=' ')
    t0 = time.perf_counter()        
    overpasses = find_overpasses(t, location, [sat], sun, min_elevation)
    tf = time.perf_counter() - t0
    if verbose: 
        print_fn(f'{tf:0.3f} sec')
    return overpasses
    


