import pickle
import datetime
from typing import List

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
from astropy import units as u
from astropy.coordinates import SkyCoord, TEME, CartesianRepresentation, ITRS, EarthLocation, AltAz
from astropy.time import Time

from .rotations.rotations import site_ECEF
from .rotations.transform import ecef2eci, ecef2sez, teme2ecef
from .rotations.polar import eop
from .solar import sun_pos, is_sat_illuminated
from .topocentric import razel, site_sat_rotations
from .propagate import propagate_satellite
from .timefn import jday2datetime, julian_date, jd2jc
from .schemas import Point, Overpass, Satellite, Location
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
    t = sat.time
    assert idx0 >= 0
    assert idxf <= len(t.jd)
    
    # Check if site is in daylight, compute dot product of sun position.
    # TODO: compute ECI vector for site position
    if np.dot(sun.rECI, rsiteECI) > 0:
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


# def get_overpasses(el, azm, rng, jdt_ary, rSEZ, min_elevation=10, sat_id=None):
#     el0 = el[:-1] - min_elevation
#     el1 = el[1:] - min_elevation
#     el_change_sign = (el0*el1 < 0)   
#     start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]  # Find the start of an overpass
#     end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]    # Find the end of an overpass
#     num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
#     if start_idx.size < end_idx.size:
#         end_idx = end_idx[1:]
#     overpasses = [None] * num_overpasses
#     for j in range(num_overpasses):
#         # Store indecies of overpasses in a list
#         idx0 = start_idx[j]
#         idxf = end_idx[j]
#         overpass_idx = np.arange(idx0, idxf+1, dtype=int)
#         idxmax = np.argmax(el[overpass_idx])
#         start_pt = Point(
#             datetime=jday2datetime(jdt_ary[idx0]),
#             azimuth=azm[idx0],
#             elevation=el[idx0],
#             range=rng[idx0]
#         )
#         max_pt = Point(
#             datetime=jday2datetime(jdt_ary[idx0 + idxmax]),
#             azimuth=azm[idx0 + idxmax],
#             elevation=el[idx0 + idxmax],
#             range=rng[idx0 + idxmax]
#         )
#         end_pt = Point(
#             datetime=jday2datetime(jdt_ary[idxf]),
#             azimuth=azm[idxf],
#             elevation=el[idxf],
#             range=rng[idxf]
#         )
#         if sat_id is not None:
#             overpass = Overpass(
#                 satellite_id=sat_id,
#                 start_pt=start_pt,
#                 max_pt=max_pt,
#                 end_pt=end_pt
#             )
#         else:
#             overpass = Overpass(
#                 start_pt=start_pt,
#                 max_pt=max_pt,
#                 end_pt=end_pt
#             )
#         overpasses[j] = overpass
#     return overpasses


# def predict_passes(lat, lon, h, rsatECEF, rsatECI, jdt, rsun=None, min_elevation=None):
#     rsiteECEF = site_ECEF(lat, lon, h)
#     rsiteECI = ecef2eci(rsiteECEF, jdt)
#     rho = site_sat_rotations(rsiteECEF, rsatECEF)
#     rSEZ = ecef2sez(rho, lat, lon)
#     rng, az, el = razel(rSEZ)
#     overpasses = get_overpasses(el, az, rng, jdt, rSEZ, rsiteECI=None, rsatECI=None, min_elevation=min_elevation)
#     return overpasses


def find_overpasses(location: Location, sats: List[SpaceObject], times: Time, sun: List[SpaceObject], min_elevation: float = 10):

    if len(sats) > 1:
        store_sat_id = True
    else:
        store_sat_id = False

    rsiteECEF = site_ECEF(location.lat, location.lon, location.h)
    overpasses = []
    
    for sat in sats:
        rho = RhoVector()
        rho.time = times
        rho.rECEF = site_sat_rotations(rsiteECEF, sat.rECEF)
        rho.rSEZ = ecef2sez(rho.rECEF, location.lat, location.lon)
        rng, az, el = razel(rho.rSEZ)
        rho.rng = rng
        rho.az = az
        rho.el = el
        # rsiteECI = ecef2eci(rsiteECEF, jdt_utc1)
        
        # Find Overpasses
        el0 = rho.el[:-1] - min_elevation
        el1 = rho.el[1:] - min_elevation
        el_change_sign = (el0*el1 < 0)   
        start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]  # Find the start of an overpass
        end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]    # Find the end of an overpass
        num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
        if start_idx.size < end_idx.size:
            end_idx = end_idx[1:]
        sat_overpasses = [None] * num_overpasses
        for j in range(num_overpasses):
            # Store indecies of overpasses in a list
            idx0 = start_idx[j]
            idxf = end_idx[j]
            overpass_idx = np.arange(idx0, idxf+1, dtype=int)
            idxmax = np.argmax(el[overpass_idx])
            start_pt = Point(
                datetime=jday2datetime(rho.time.jd[idx0]),
                azimuth=rho.az[idx0],
                elevation=rho.el[idx0],
                range=rho.rng[idx0]
            )
            max_pt = Point(
                datetime=jday2datetime(rho.time.jd[idx0 + idxmax]),
                azimuth=rho.az[idx0 + idxmax],
                elevation=rho.el[idx0 + idxmax],
                range=rho.rng[idx0 + idxmax]
            )
            end_pt = Point(
                datetime=jday2datetime(rho.time.jd[idxf]),
                azimuth=rho.az[idxf],
                elevation=rho.el[idxf],
                range=rho.rng[idxf]
            )
            if store_sat_id:
                overpass = Overpass(
                    satellite_id=satellite.id,
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
            sat_overpasses[j] = overpass

            # if verbose:
            #     print('Determine visibility')

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
        dt_start = datetime.datetime.now()
    if dt_end is None:
        dt_end = dt_start + datetime.timedelta(days=14)
    if tle is None:
        tle = get_TLE(satellite)

    jdt0 = julian_date(dt_start)
    jdtf = julian_date(dt_end)
    total_days = (dt_start-dt_end).total_seconds()/60
    dt_days = dt_seconds/(24*60*60.0)
    jd_array = np.arange(jdt0, jdtf, dt_days, dtype=float)
    t = Time(jd_array, format='jd')
    # t.tt = jd2jc(t.jd)
    # dUTC1, xp, yp = eop(t.jd)
    # t.jd_utc1 = t.jd + dUTC1
    # t.tt_utc1 = jd2jc(t.jd_utc1)

    sat = SpaceObject()
    sat.time = t

    if verbose:
        print(f"begin propagation from {dt_start.isoformat()} to {dt_end.isoformat()}...")
    rTEME, _ = propagate_satellite(tle.tle1, tle.tle2, t.jd)

    # Use the TEME reference frame from astropy
    rTEME_astropy = CartesianRepresentation(rTEME * u.km)
    teme = TEME(rTEME_astropy, obstime=t)
    sat.rECEF_astropy = teme.transform_to(ITRS(obstime=t))
    sat.subpoint = sat.rECEF_astropy.earth_location
    
    if verbose:
        print(f"rotate satellite position from TEME to ECEF...")
    sat.rECEF = teme2ecef(rTEME, t.jd + t.delta_ut1_utc)
    sat.rECI = rTEME.view()

    # Compute sun-satellite quantities
    if verbose:
        print(f"Compute sun-satellite quantities...")
    sun = Sun()
    sun.time = t
    sun.rECI = sun_pos(sun.time.jd)  # to do: use cached value
    sat.illuminated = is_sat_illuminated(rTEME.copy(), sun.rECI)
        
    if verbose:
        print('begin prediction...')
    overpasses = find_overpasses(location, [sat], t, sun, min_elevation)
    return overpasses
    


