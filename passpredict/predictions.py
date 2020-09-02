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
    date_start: datetime = None,
    date_end: datetime = None,
    dt_seconds: int = 1,
    min_elevation: float = 10.0,
    cache: 'Cache' = None,
    verbose: bool = False, 
    store_sat_id: bool = False, 
    print_fn: Callable = print, 
    visible_only: bool = False
):
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

    jd = julian_date_array_from_date(date_start, date_end, dt_seconds)

    if cache is None:
        t = Time(jd, format='jd')
        sun = compute_sun_data(t)
        tle = get_TLE(satellite.id)
        sat = compute_satellite_data(tle, t, sun)

    else:

        with cache:

            time_key = 'time:' + date_start.strftime('%Y%m%d') + date_end.strftime('%Y%m%d') + str(dt_seconds)

            sun_key = 'sun:' + time_key
            sun = cache.get(sun_key)
            if sun is None:
                if verbose:
                    print_fn("Compute sun position... ")
                t = Time(jd, format='jd')
                sun = compute_sun_data(t)
                cache.set(sun_key, sun, ttl=86400)

            sat_key = 'sat:' + str(satellite.id) + time_key
            sat = cache.get(sat_key)
            if sat is None:
                if verbose:
                    print_fn(f"Propagate satellite from {date_start.isoformat()} to {date_end.isoformat()}... ")
                tle = get_TLE(satellite.id)
                t = Time(jd, format='jd')
                sat = compute_satellite_data(tle, t, sun)
                cache.set(sat_key, sat, ttl=86400)

    if verbose:
        print_fn('Predict overpasses... ')
    overpasses = find_overpasses(jd, location, [sat], sun, min_elevation, visible_only=visible_only)
    return overpasses



