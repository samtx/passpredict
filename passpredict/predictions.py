import pickle
from datetime import datetime, timedelta
from typing import List

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
from astropy import units as u
from astropy.coordinates import SkyCoord, TEME, CartesianRepresentation, ITRS, EarthLocation, AltAz, get_sun
from astropy.time import Time

from .rotations.rotations import site_ECEF
from .rotations.transform import ecef2eci, ecef2sez, teme2ecef
from .rotations.polar import eop
from .solar import sun_pos, is_sat_illuminated
from .topocentric import razel, site_sat_rotations
from .propagate import propagate_satellite
from .timefn import jday2datetime, julian_date, jd2jc
from .schemas import Point, Overpass, Satellite, Location, Tle
from .models import SatelliteRV, SpaceObject, Sun, RhoVector
from .utils import get_TLE


# Reference: https://docs.astropy.org/en/latest/coordinates/satellites.html

def vector_angle(r1, r2):
    """Compute the angle between two vectors

        r1 and r2 : (3, n), n is the number of observations
    """
    numerator = np.einsum("ij,ij->j", r1, r2)  # vectorized dot product
    denominator = norm(r1, axis=0) * norm(r2, axis=0)
    out = np.arccos(numerator / denominator) * RAD2DEG
    return out


def determine_visibilty(idx0: int, idxf: int, sat: SpaceObject, loc, rho: RhoVector, sun: Sun):
    """
    Determine satellite visibility
    """
    assert sat.time == sun.time == rho.time
    assert idx0 >= 0
    assert idxf <= len(t.jd)
    t = sat.time
    
    # Check if site is in daylight, compute dot product of sun position.
    if np.dot(sun.rECEF, rsiteECI) > 0:
        # site is in daylight
        return 1
    else:
        # If nighttime, check if satellite is illuminated or in shadow
        if is_sat_illuminated(rsatECI_i, rsun):
            # satellite is illuminated
            return 3
        else:
            # satellite is in shadow
            return 2


def satellite_visible(rsatECI, rsiteECI, rho, jdt):
    """Determine visibility of satellite from Earth"""
    visible = np.zeros(jdt.size)
    # First, determine if satellite is above the horizon
    # find indecies where rho[2] > 0 --> idx
    vis_idx = np.nonzero(rho[2] > 0)[0]
    # Loop over indecies for times that satellite is over horizon
    for i in range(len(vis_idx)):
        idx = vis_idx[i]
        jdt_i = jdt[idx]
        rsatECI_i = rsatECI[:, idx]
        rsiteECI_i = rsiteECI[:, idx]
        rho_i = rho[:, idx]
        # Check if site is in daylight, compute dot product of sun position.
        rsun = sun_pos(jdt_i)
        if len(rsun.shape) > 1:
            rsun = rsun.flatten()
        # TODO: compute ECI vector for site position
        if np.dot(rsun, rsiteECI) > 0:
            # site is in daylight
            visible[idx] = 1
        else:
            # If nighttime, check if satellite is illuminated or in shadow
            if is_sat_illuminated(rsatECI_i, rsun):
                # satellite is illuminated
                visible[idx] = 3
            else:
                # satellite is in shadow
                visible[idx] = 2
    return visible



def find_overpasses(location: Location, sats: List[SpaceObject], times: Time, sun: List[SpaceObject], min_elevation: float = 10) -> List[Overpass]:
    """
    Real-time computation for finding satellite overpasses of a topographic location.
    Can support multiple satellites over a single location
    """
    store_sat_id = True if len(sats) > 0 else False
    overpasses = []
    for sat in sats:
        rho = RhoVector(sat, location)
        sat_overpasses = rho.find_overpasses(min_elevation, store_sat_id)
        overpasses += sat_overpasses

    return overpasses

def compute_time_array(dt_start: datetime, dt_end: datetime, dt_seconds: float) -> Time:
    """
    Create astropy Time object
    """
    jdt0 = julian_date(dt_start)
    jdtf = julian_date(dt_end)
    total_days = (dt_start-dt_end).total_seconds()/60
    dt_days = dt_seconds/(24*60*60.0)
    jd_array = np.arange(jdt0, jdtf, dt_days, dtype=float)
    return Time(jd_array, format='jd')


def compute_satellite_data(tle: Tle, t: Time) -> SpaceObject:
    """
    Compute satellite data for Time
    """
    sat = SpaceObject()
    sat.time = t
    r, _ = propagate_satellite(tle.tle1, tle.tle2, t.jd)
    # Use the TEME reference frame from astropy
    teme = TEME(CartesianRepresentation(r * u.km), obstime=t)
    ecef = teme.transform_to(ITRS(obstime=t))
    sat.rECEF = ecef.data.xyz.value  # extract numpy array from astropy object
    sat.subpoint = ecef.earth_location
    sat.latitude = sat.subpoint.lat.value
    sat.longitude = sat.subpoint.lon.value
    # sat.meta.id = tle.satellite.id
    return sat
    

def compute_sun_data(t: Time) -> Sun:
    """
    Compute sun position data

    Compute for each minute, then interpolate for each second
    """
    t_tmp = t[::60]
    sun_tmp = get_sun(t_tmp)  # get astropy coordinates for sun in GCRS
    sun_tmp = sun_tmp.transform_to(ITRS(obstime=t_tmp))  # transform to ECEF frame
    sun_tmp = sun_tmp.data.xyz.to('km').value
    sun_data = np.empty((3, t.size))
    for i in range(3):
        sun_data[i] = np.interp(t.jd, t_tmp.jd, sun_tmp[i])
    sun = Sun()
    sun.time = t
    sun.rECEF = sun_data
    return sun


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
        print(f"Begin propagation from {dt_start.isoformat()} to {dt_end.isoformat()}...")
    sat = compute_satellite_data(tle, t)

    if verbose:
        print("Compute sun position...")
    sun = compute_sun_data(t)

    # Compute sun-satellite quantities
    if verbose:
        print("Compute sun-satellite illumination...")
    sat.illuminated = is_sat_illuminated(sat.rECEF, sun.rECEF)

        
    if verbose:
        print('begin prediction...')
    overpasses = find_overpasses(location, [sat], t, sun, min_elevation)
    return overpasses
    


