import pickle
import datetime

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm

from .rotations import site_sat_rotations
from .solar import sun_pos, is_sat_illuminated
from .topocentric import razel
from .propagate import propagate
from .timefn import jday2datetime
from .schemas import Point, Overpass, SatelliteRV, Satellite
from .utils import get_TLE

def vector_angle(r1, r2):
    """Compute the angle between two vectors

        r1 and r2 : (3, n), n is the number of observations
    """
    numerator = np.einsum("ij,ij->j", r1, r2)  # vectorized dot product
    denominator = norm(r1, axis=0) * norm(r2, axis=0)
    out = np.arccos(numerator / denominator) * RAD2DEG
    return out


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


def get_overpasses(el, azm, rng, jdt_ary, rSEZ, rsiteECI=None, rsatECI=None, min_elevation=10, sat_id=None):
    el0 = el[:-1] - min_elevation
    el1 = el[1:] - min_elevation
    el_change_sign = (el0*el1 < 0)   
    start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]  # Find the start of an overpass
    end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]    # Find the end of an overpass
    num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
    if start_idx.size < end_idx.size:
        end_idx = end_idx[1:]
    # overpasses = np.empty(num_overpasses, dtype=object)
    overpasses = [None] * num_overpasses
    for j in range(num_overpasses):
        # Store indecies of overpasses in a list
        idx0 = start_idx[j]
        idxf = end_idx[j]
        overpass_idx = np.arange(idx0, idxf+1, dtype=int)
        idxmax = np.argmax(el[overpass_idx])
        start_pt = Point(
            datetime=jday2datetime(jdt_ary[idx0]),
            azimuth=azm[idx0],
            elevation=el[idx0],
            range=rng[idx0]
        )
        max_pt = Point(
            datetime=jday2datetime(jdt_ary[idx0 + idxmax]),
            azimuth=azm[idx0 + idxmax],
            elevation=el[idx0 + idxmax],
            range=rng[idx0 + idxmax]
        )
        end_pt = Point(
            datetime=jday2datetime(jdt_ary[idxf]),
            azimuth=azm[idxf],
            elevation=el[idxf],
            range=rng[idxf]
        )
        # sat_vis = satellite_visible(rsatECI, rsiteECI, rSEZ, jdt)
        if sat_id is not None:
            overpass = Overpass(
                satellite_id=sat_id,
                start_pt=start_pt,
                max_pt=max_pt,
                end_pt=end_pt
            )
        else:
            overpass = Overpass(
                start_pt=start_pt,
                max_pt=max_pt,
                end_pt=end_pt
            )
        overpasses[j] = overpass
    return overpasses


def predict_passes(lat, lon, h, rsatECEF, rsatECI, jdt, rsun=None, min_elevation=None):
    rSEZ = site_sat_rotations(lat, lon, h, rsatECEF)
    # rsiteECI = site2eci(lat, lon, h, jdt)
    rng, az, el = razel(rSEZ)
    # plot_elevation(np.arange(el.size), el)
    overpasses = get_overpasses(el, az, rng, jdt, rSEZ, rsiteECI=None, rsatECI=None, min_elevation=min_elevation)
    return overpasses


def predict(location, satellite, dt_start=None, dt_end=None, dt_seconds=1, min_elevation=None, tle=None, reload=True):
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
        dt_start = datetime.datetime.now()
    if dt_end is None:
        dt_end = dt_start + datetime.timedelta(days=14)
    if tle is None:
        tle = get_TLE(satellite)
    print(f"begin propagation from {dt_start.isoformat()} to {dt_end.isoformat()}")
    _reload = reload
    if _reload:
        satellite_rv = propagate(tle.tle1, tle.tle2, dt_start, dt_end, dt_seconds)
        satellite_rv.satellite = satellite
        satellite_rv.tle = tle
        # Compute sun-satellite quantities
        jdt = satellite_rv.julian_date
        rsunECI = sun_pos(jdt)
        satellite_rv.visible = is_sat_illuminated(satellite_rv.rECI, rsunECI)
        
        with open(f'satellite_{satellite.id:d}.pkl', 'wb') as f:
            pickle.dump(satellite_rv, f)
    else:
        with open(f'satellite_{satellite.id:d}.pkl', 'rb') as f:
            satellite_rv = pickle.load(f)
    print('begin prediction...')
    # set minimum elevation parameter: min_elevation = 10 degrees
    overpasses = predict_passes(
        location.lat, location.lon, location.h,
        satellite_rv.rECEF, satellite_rv.rECI, satellite_rv.julian_date,
        min_elevation=min_elevation)
    return overpasses


