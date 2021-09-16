from __future__ import annotations
from datetime import timedelta, date
from typing import List, NamedTuple
import math

from numpy import ndarray
import numpy as np
from fastapi import HTTPException

from app.schemas import Overpass, Location, Satellite, OverpassResult, Point, PassType, Tle
from app.models import SatPredictData
from app.settings import MAX_DAYS, DT_SECONDS
from app.utils import get_visible_satellites
from app.tle import get_most_recent_tle
from app.astrodynamics import _overpass
from app.astrodynamics._solar import sun_pos_ecef
from app.astrodynamics.propagate import compute_satellite_data
from app.astrodynamics.timefn import julian_date_array_from_date, jday2datetime
from app.astrodynamics.constants import RAD2DEG, DAY_S
from app.astrodynamics._rotations import ecef2sez
from app.astrodynamics.topocentric import site_ECEF


VISIBLE_SATS = get_visible_satellites()

# Make sure that dt_seconds evenly divides into number of seconds per day
assert DAY_S % DT_SECONDS == 0
NUM_TIMESTEPS_PER_DAY = int(DAY_S / DT_SECONDS)
TIME_KEY_SUFFIX = ':' + str(MAX_DAYS) + ':' + str(DT_SECONDS)


class RhoVector:
    def __init__(self, jd: np.ndarray, object_rECEF: np.ndarray, location: Location):
        self.jd = jd
        self.location = location
        self.site_ECEF = site_ECEF(location.lat, location.lon, location.h)[np.newaxis, :]
        self.rECEF = object_rECEF - self.site_ECEF
        self.rSEZ = ecef2sez(self.rECEF, location.lat, location.lon)
        self.rng, self.el = self._compute_range_and_elevation()

    def _compute_range_and_elevation(self):
        return _overpass.compute_range_and_elevation(self.rSEZ)

    def _azimuth_from_idx(self, idx):
        rS, rE = self.rSEZ[idx, [0,1]]
        tmp = np.arctan2(rS, rE)
        az = (tmp + math.pi * 0.5) * RAD2DEG
        if rS < 0 and rE < 0:
            az %= 360
        return az

    def point(self, idx):
        dt = jday2datetime(self.jd[idx])
        p = Point.construct(
            datetime=dt,
            timestamp=dt.timestamp(),
            azimuth=round(self._azimuth_from_idx(idx), 2),
            elevation=round(self.el[idx], 2),
            range=round(self.rng[idx], 3)
        )
        return p

    def __getitem__(self, slc):
        return RhoVector(
            jd=self.jd[slc],
            object_rECEF=self.rECEF[slc, :],
            location=self.location,
        )


class SatelliteVisibility(NamedTuple):
    passtype: PassType
    vis_start_pt: Point = None
    vis_end_pt: Point = None
    brightness: float = None


def compute_satellite_visibility(
    sun_rho: RhoVector,
    sat_rho: RhoVector,
    sunset_el: float,
    sat: SatPredictData
) -> SatelliteVisibility:
    """
    Function to compute whether the satellite is visible to the observer
    """
    assert sun_rho.location == sat_rho.location
    assert sun_rho.rECEF.shape == sat_rho.rECEF.shape
    assert sun_rho.rSEZ.shape == sat_rho.rSEZ.shape
    np.testing.assert_allclose(sun_rho.jd, sat_rho.jd)
    site_in_sunset = sun_rho.el - sunset_el < 0
    site_in_sunset_idx = np.nonzero(site_in_sunset)[0]
    if site_in_sunset_idx.size == 0:
        passtype = PassType.daylight # site is always sunlit, so overpass is in daylight
    else:
        # get satellite illumination values for this overpass
        sat_illuminated = sat.sun_sat_dist > 0
        sat_visible = (sat_illuminated * site_in_sunset)
        if np.any(sat_visible):
            passtype = PassType.visible  # site in night, sat is illuminated
            sat_visible_idx = np.nonzero(sat_visible)[0]
            sat_visible_start_idx = sat_visible_idx.min()
            sat_visible_end_idx = sat_visible_idx.max()
            sat_rho_vis = sat_rho[sat_visible_start_idx:sat_visible_end_idx + 1]
            sun_rho_vis = sun_rho[sat_visible_start_idx:sat_visible_end_idx + 1]
            vis_start_pt = sat_rho_vis.point(0)
            vis_end_pt = sat_rho_vis.point(-1)
            brightness_idx = np.argmax(sat_rho_vis.el)
            sat_site_sun_angle = math.acos(
                np.dot(sat_rho_vis.rSEZ[brightness_idx,:], sun_rho_vis.rSEZ[brightness_idx,:]) \
                / (sat_rho_vis.rng[brightness_idx] * sun_rho_vis.rng[brightness_idx])
            )
            beta = math.pi - sat_site_sun_angle  # phase angle: site -- sat -- sun angle
            brightness = sat.intrinsic_mag - 15 + 5*math.log10(sat_rho_vis.rng[brightness_idx]) \
                - 2.5*math.log10(math.sin(beta) + (math.pi - beta)*math.cos(beta))
        else:
            passtype = PassType.unlit  # nighttime, not illuminated (radio night)
    if passtype != PassType.visible:
        vis_start_pt = vis_end_pt = brightness = None
    satellite_visitbility = SatelliteVisibility(
        passtype=passtype, vis_start_pt=vis_start_pt, vis_end_pt=vis_end_pt, brightness=brightness
    )
    return satellite_visitbility


def compute_single_satellite_overpasses(
    sat,
    *,
    jd=None,
    location=None,
    sun_rECEF=None,
    min_elevation=10.0,
    visible_only=False,
    store_sat_id=True,
    sunset_el=-8.0
) -> List[Overpass]:
    rho = RhoVector(jd, sat.rECEF, location)
    start_idx, end_idx = _start_end_index(rho.el - min_elevation)
    num_overpasses = min(start_idx.size, end_idx.size)
    if start_idx.size < end_idx.size:
        end_idx = end_idx[1:]
    overpasses = [None]*num_overpasses if not visible_only else []
    for j in range(num_overpasses):
        # Store indecies of overpasses in a list
        idx0 = start_idx[j]
        idxf = end_idx[j]
        idxmax = np.argmax(rho.el[idx0:idxf+1])
        start_pt = rho.point(idx0)
        max_pt = rho.point(idx0 + idxmax)
        end_pt = rho.point(idxf)
        # Find visible start and end times
        if sun_rECEF is not None:
            sun_rho = RhoVector(
                jd=jd[idx0:idxf+1],
                object_rECEF=sun_rECEF[idx0:idxf+1],
                location=location
            )
            satellite_visibility = compute_satellite_visibility(
                sun_rho=sun_rho, sat_rho=rho[idx0:idxf+1],
                sunset_el=sunset_el, sat=sat[idx0:idxf+1]
            )
            passtype = satellite_visibility.passtype
        else:
            passtype = None
        if visible_only and passtype != PassType.visible:
            continue
        overpass_dict = {
            'start_pt': start_pt,
            'max_pt': max_pt,
            'end_pt': end_pt,
        }
        if store_sat_id:
            overpass_dict['satellite_id'] = sat.id
        if (passtype is not None) and (passtype == PassType.visible):
            overpass_dict['vis_start_pt'] = satellite_visibility.vis_start_pt
            overpass_dict['vis_end_pt'] = satellite_visibility.vis_end_pt
            overpass_dict['brightness'] = round(satellite_visibility.brightness, 2)
        overpass_dict['type'] = passtype
        overpass = Overpass.construct(**overpass_dict)
        if not visible_only:
            overpasses[j] = overpass
        else:
            overpasses.append(overpass)
    return overpasses


def predict_single_satellite_overpasses(
    satid: int,
    location: Location,
    *,
    date_start: date = None,
    days: int = 10,
    min_elevation: float = 10.0,
    visible_only: bool = False,
    db: Connection = None,  # database connection
    cache: Redis.client = None,   # cache connection
) -> OverpassResult:
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
        date_start = date.today()

    # Get TLE data
    tle_key = f"tle:{satid}:{date_start.strftime('%Y%m%d')}"
    result = cache.hgetall(tle_key)
    if result:
        tle1 = result[b'tle1'].decode()
        tle2 = result[b'tle2'].decode()
        tle = Tle.from_string(tle1=tle1, tle2=tle2)
        tle_is_new = False
    else:
        tle = get_most_recent_tle(db, satid, raise_404=True)
        data = {
            'tle1': tle.tle1,
            'tle2': tle.tle2
        }
        cache.hset(tle_key, mapping=data)
        cache.expire(tle_key, 28800)  # save tle in cache for 8 hours
        tle_is_new = True

    date_end = date_start + timedelta(days=days)
    jd = julian_date_array_from_date(date_start, date_end, DT_SECONDS)
    n = jd.size
    time_key = f"time:{date_start.strftime('%Y%m%d')}:{date_end.strftime('%Y%m%d')}:{DT_SECONDS}"

    # Get Sun position data
    sun_key = f'sun:{time_key}'
    result = cache.get(sun_key)
    if result:
        sun_rECEF = np.frombuffer(result, dtype=np.float64).reshape((n, 3))
    else:
        sun_rECEF = sun_pos_ecef(jd)
        cache.set(sun_key, sun_rECEF.tobytes(), ex=86400)

    sat_key = f'sat:{satid}:{time_key}'
    if (not tle_is_new) and (result := cache.hgetall(sat_key)):
        sat_rECEF = np.frombuffer(result[b'rECEF'], dtype=np.float64).reshape((n, 3))
        sun_sat_dist = np.frombuffer(result[b'sun_sat_dist'], dtype=np.float64)
        sat = SatPredictData(id=satid, rECEF=sat_rECEF, sun_sat_dist=sun_sat_dist)
    else:
        # if we got a new TLE, then always recompute
        sat = compute_satellite_data(tle, jd, sun_rECEF)
        data = {
            'rECEF': sat.rECEF.tobytes(),
            'sun_sat_dist': sat.sun_sat_dist.tobytes()
        }
        cache.hset(sat_key, mapping=data)
        cache.expire(sat_key, 86400)

    overpasses = compute_single_satellite_overpasses(
        sat,
        jd=jd,
        location=location,
        sun_rECEF=sun_rECEF,
        min_elevation=min_elevation,
        visible_only=visible_only,
        store_sat_id=False
    )
    overpass_result = OverpassResult(
        location=location,
        satid=satid,
        overpasses=overpasses
    )
    return overpass_result


def predict_all_visible_satellite_overpasses(
    location: Location,
    *,
    date_start: date = None,
    min_elevation: float = 10.0,
    db: Connection = None,
    cache: Redis.Client = None
):
    if date_start is None:
        date_start = date.today()
    date_end = date_start + timedelta(days=1)
    jd = julian_date_array_from_date(date_start, date_end, DT_SECONDS)
    sun_rECEF = sun_pos_ecef(jd)

    # Query TLEs for visible satellites
    try:
        conn = engine.connect()
        s = select([tledb.c.satellite_id, tledb.c.tle1, tledb.c.tle2]) \
            .where(tledb.c.satellite_id.in_(VISIBLE_SATS))
        res = conn.execute(s)

    except:
        print("Error quering TLEs")
    finally:
        conn.close()

    overpasses = []
    for satid in VISIBLE_SATS:
        tle = get_most_recent_tle(db, satid, raise_404=False)
        if not tle:
            # no TLE data for satellite
            continue
        sat = compute_satellite_data(tle, jd, sun_rECEF)
        sat_overpasses = compute_single_satellite_overpasses(
            sat,
            jd=jd,
            location=location,
            sun_rECEF=sun_rECEF,
            min_elevation=min_elevation,
            visible_only=True,
            store_sat_id=True
        )
        overpasses += sat_overpasses
    overpass_result = OverpassResult(
        location=location,
        overpasses=overpasses
    )
    return overpass_result


def _start_end_index(x):
    """
    Finds the start and end indecies when a 1D array crosses zero
    """
    x0 = x[:-1]
    x1 = x[1:]
    x_change_sign = (x0*x1 < 0)
    start_idx = np.nonzero(x_change_sign & (x0 < x1))[0]  # Find the start of an overpass
    end_idx = np.nonzero(x_change_sign & (x0 > x1))[0]    # Find the end of an overpass
    return start_idx, end_idx


def set_sun_cache(
    sun_key: str,
    rECEF: np.ndarray,
    pipe: Pipeline,
    ttl: int = 86400
):
    """
    Add sun hashdata to redis pipeline
    """
    rECEF_value = {
        'n': rECEF.shape[0],
        'dtype': str(rECEF.dtype),
        'array': rECEF.tobytes()
    }
    pipe.hset(sun_key, mapping=rECEF_value)
    pipe.expire(sun_key, ttl)
    return pipe


def get_sun_cache(sun: bytes):
    """
    Get numpy array from Redis bytes
    """
    n = int(sun[b'n'])
    dtype = sun[b'dtype'].decode()
    arr = np.frombuffer(sun[b'array'], dtype=dtype).reshape((n, 3))
    return arr


def set_sat_cache(sat_key: str, sat: SatPredictData, pipe: Pipeline, ttl: int=86400):
    """
    Set satellite hashtable cache data to pipeline
    """
    sat_rECEF = sat.rECEF
    data = {
        'n': sat_rECEF.shape[0],
        'rECEF': sat_rECEF.tobytes(),
        'sun_sat_dist': sat.sun_sat_dist.tobytes()
    }
    pipe.hset(sat_key, mapping=data)
    pipe.expire(sat_key, ttl)
    return pipe


def get_sat_cache(sat: bytes, satid: int):
    n = int(sat[b'n'])
    rECEF = np.frombuffer(sat[b'rECEF'], dtype=np.float64).reshape((n, 3))
    sun_sat_dist = np.frombuffer(sat[b'sun_sat_dist'], dtype=np.float64)
    return SatPredictData(id=satid, rECEF=rECEF, sun_sat_dist=sun_sat_dist)
