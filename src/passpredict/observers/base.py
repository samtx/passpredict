from __future__ import annotations
import datetime
from math import pi, log10, sin, cos, acos, degrees, radians, floor
import typing
from functools import lru_cache, cached_property
from collections import defaultdict
from dataclasses import dataclass, asdict
from enum import Enum
from abc import abstractmethod

import numpy as np
from numpy.linalg import norm
from orbit_predictor.predictors.pass_iterators import LocationPredictor

from .functions import julian_date_sum
from ..time import julian_date_from_datetime
from ..constants import R_EARTH
from .. import _rotations
from ..utils import get_pass_detail_datetime_metadata
from ..exceptions import NotReachable
from ..solar import sun_pos

if typing.TYPE_CHECKING:
    from ..satellites import LLH
    from ..locations import Location


class RangeAzEl(typing.NamedTuple):
    range: float  # km
    az: float     # deg
    el: float     # deg


@dataclass(frozen=True)
class PassPoint:
    dt: datetime                # datetime UTC
    range: float                # km
    azimuth: float              # deg
    elevation: float            # deg
    brightness: float = None    # magnitude brightness
    type: Visibility = None       # visibility status

    @cached_property
    def direction(self) -> str:
        ''' Return ordinal direction from azimuth degree '''
        azm = self.azimuth % 360
        mod = 360/16. # number of degrees per coordinate heading
        start = 0 - mod/2
        n = int(floor((azm-start)/mod))
        coordinates = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','N']
        return coordinates[n]


class PassType(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'


class Visibility(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'
    not_visible = 'not visible'


class BasicPassInfo:
    """
    Holds basic pass information:
        aos datetime
        los datetime
        tca datetime
        max elevation
        duration (property)
        valid (property)
    """
    def __init__(
        self,
        aos_dt: datetime.datetime,
        tca_dt: datetime.datetime,
        los_dt: datetime.datetime,
        max_elevation: float,
        type_: PassType = None,
        vis_begin_dt: datetime.datetime = None,
        vis_end_dt: datetime.datetime = None,
        vis_tca_dt: datetime.datetime = None,
    ):
        self.aos_dt = aos_dt if aos_dt is not None else None
        self.tca_dt = tca_dt
        self.los_dt = los_dt if los_dt is not None else None
        self.max_elevation = max_elevation
        self.type = type_
        self.vis_begin_dt = vis_begin_dt
        self.vis_end_dt = vis_end_dt
        self.vis_tca_dt = vis_tca_dt

    @property
    def aos(self):
        """ Backwards compatibilty from orbit_predictor """
        return self.aos_dt

    @property
    def tca(self):
        """ Backwards compatibilty from orbit_predictor """
        return self.tca_dt

    @property
    def los(self):
        """ Backwards compatibilty from orbit_predictor """
        return self.los_dt

    @property
    def valid(self):
        if (self.aos is None) or (self.los is None):
            return False
        return (self.max_elevation > 0)

    @cached_property
    def max_elevation_deg(self):
        return degrees(self.max_elevation)

    @cached_property
    def duration(self) -> datetime.timedelta:
        return self.los - self.aos


@dataclass
class PredictedPass:
    satid: int
    location: Location
    aos: PassPoint
    tca: PassPoint
    los: PassPoint
    type: PassType = None
    azimuth: typing.Sequence[float] = None
    elevation: typing.Sequence[float] = None
    range: typing.Sequence[float] = None
    datetime: typing.Sequence[datetime.datetime] = None
    vis_begin: PassPoint = None
    vis_end: PassPoint = None
    vis_tca: PassPoint = None
    brightness: float = None

    @cached_property
    def midpoint(self):
        """ Return datetime in UTC of midpoint of pass """
        midpt = self.aos.dt + (self.los.dt - self.aos.dt) / 2
        return midpt

    @cached_property
    def duration(self) -> float:
        """ Return pass duration in seconds """
        return (self.los.dt - self.aos.dt).total_seconds()

    def __repr__(self):
        return f"<PredictedPass {self.satid} over {repr(self.location)} on {self.aos.dt}"

    def dict(self):
        """
        Serialize into dictionary
        """
        data = {
            'satid': self.satid,
            'location': {
                'name': self.location.name,
                'lat': self.location.latitude_deg,
                'lon': self.location.longitude_deg,
                'h': self.location.elevation_m,
            },
            'aos': asdict(self.aos),
            'tca': asdict(self.tca),
            'los': asdict(self.los),
            'type': self.type
        }
        return data


class ObserverBase(LocationPredictor):

    def __init__(
        self,
        location,
        satellite,
        max_elevation_gt: float = 0,
        aos_at_dg: float = 0,
        tolerance_s: float = 1,
        sunrise_dg: float = -6,
    ):
        """
        Initialize Observer but also compute radians for geodetic coordinates
        """
        self.location = location
        self.satellite = satellite
        self.max_elevation_gt = radians(max([max_elevation_gt, aos_at_dg]))
        self.aos_at = radians(aos_at_dg)
        self.aos_at_dg = aos_at_dg
        if tolerance_s <= 0:
            raise Exception("Tolerance must be > 0")
        self.tolerance_s = tolerance_s
        self.tolerance = datetime.timedelta(seconds=tolerance_s)
        self.jd_tol = self.tolerance_s / 86400
        self.sunrise_dg = sunrise_dg

    @property
    def predictor(self):
        """ For backwards compatibility with orbit-predictor location predictor """
        return self.satellite

    @abstractmethod
    def iter_passes(self, start_date, limit_date=None):
        raise NotImplementedError

    def _is_pass_valid(self, pass_, visible_only=False):
        if (pass_.aos is None) or (pass_.los is None):
            return False
        if visible_only and pass_.type != PassType.visible:
            return False
        return (pass_.max_elevation > 0)

    @property
    def passes_over(self, *a, **kw):
        return self.iter_passes(*a, **kw)

    def get_next_pass(self,
        aos_dt: datetime.datetime,
        *,
        limit_date: datetime.datetime = None,
        visible_only: bool = False,
    ) -> PredictedPass:
        """
        Gets first overpass starting at aos_dt
        """
        pass_ = next(self.iter_passes(aos_dt, limit_date=limit_date, visible_only=visible_only), None)
        if not pass_:
            raise NotReachable('Propagation limit date exceeded')
        return pass_

    def get_next_pass_detail(
        self,
        aos_dt: datetime.datetime,
        *,
        limit_date: datetime.datetime = None,
        delta_s: float = 10,
        pad_minutes: int = 5,
    ) -> typing.Tuple[PredictedPass, LLH]:
        """
        Add details to PredictedPass
        Evaluate position and velocity properties for each delta_s seconds
        """
        pass_ = self.get_next_pass(aos_dt, limit_date=limit_date)
        start_date, n_steps, time_step = get_pass_detail_datetime_metadata(pass_, delta_s, pad_minutes=pad_minutes)
        pass_detail = self._get_overpass_detail(pass_, start_date, n_steps, time_step)
        llh = self.satellite.get_position_detail(start_date, n_steps, time_step)
        return (pass_detail, llh)

    def _get_overpass_detail(
        self,
        pass_,
        start_date: datetime.datetime,
        n_steps: int,
        time_step: datetime.timedelta,
    ):
        pass_.azimuth = np.empty(n_steps)
        pass_.elevation = np.empty(n_steps)
        pass_.range = np.empty(n_steps)
        pass_.datetime = [None] * n_steps
        dt = start_date
        for i in range(n_steps):
            rae = self.razel(dt)
            pass_.datetime[i] = dt
            pass_.azimuth[i] = rae.az
            pass_.elevation[i] = rae.el
            pass_.range[i] = rae.range
            dt += time_step
        return pass_

    def _build_predicted_pass(self, basic_pass: BasicPassInfo):
        """Returns a classic predicted pass"""
        data = defaultdict()
        data.update({
            'satid': self.satellite.satid,
            'location': self.location,
            'type': basic_pass.type,
            'aos': self.point(basic_pass.aos_dt),
            'tca': self.point(basic_pass.tca_dt),
            'los': self.point(basic_pass.los_dt),
        })
        if basic_pass.vis_begin_dt:
            data['vis_begin'] = self.point(basic_pass.vis_begin_dt)
        if basic_pass.vis_end_dt:
            data['vis_end'] = self.point(basic_pass.vis_end_dt)
        return PredictedPass(**data)

    @lru_cache(maxsize=16)
    def _elevation_at(self, when_utc):
        position = self.satellite.get_only_position(when_utc)
        return self.location.elevation_for(position)

    @lru_cache(maxsize=16)
    def _elevation_at_jd(self, jd: float) -> float:
        position = self.satellite.get_only_position_jd(jd)
        return self.location.elevation_for(position)

    def razel(self, datetime: datetime.datetime) -> RangeAzEl:
        """
        Get range, azimuth, and elevation for datetime
        """
        satellite_ecef = self.satellite.get_only_position(datetime)
        range_, az, el = _rotations.razel(
            self.location.latitude_rad, self.location.longitude_rad, self.location.recef, satellite_ecef
        )
        return RangeAzEl(range_, az, el)

    @lru_cache(maxsize=16)
    def range_at_jd(self, jd: float) -> float:
        """
        Get slant range magnitude from location to satellite [km]
        """
        satellite_ecef = self.satellite.get_only_position_jd(jd)
        range_ = _rotations.range_at(
            self.location.latitude_rad,
            self.location.longitude_rad,
            self.location.recef,
            satellite_ecef,
        )
        return range_

    def point(self, d: datetime.datetime, visibility: bool = True) -> PassPoint:
        """
        Get PassPoint with range, azimuth, and elevation data for datetime
        """
        rnazel = self.razel(d)
        vis_state = self.determine_visibility(d)
        if (vis_state == Visibility.visible) and (self.satellite.intrinsic_mag is not None):
            jd = julian_date_sum(d)
            brightness = self.brightness(jd)
        else:
            brightness = None
        pt = PassPoint(d, rnazel.range, rnazel.az, rnazel.el, type=vis_state, brightness=brightness)
        return pt

    @lru_cache(maxsize=16)
    def rho_jd(self, jd: float) -> np.ndarray:
        """
        Get topocentric ECEF vector from location to satellite
        """
        rho = np.empty(3, dtype=np.double)
        sat_recef = self.satellite.get_only_position_jd(jd)
        _rotations.ecef_to_rhosez(self.location.latitude_rad, self.location.longitude_rad, self.location.recef, sat_recef, rho)
        return rho

    def brightness(self, jd: float) -> float:
        beta = pi - self.sat_location_sun_angle(jd)
        range_ = norm(self.rho_jd(jd))
        mag = self.satellite.intrinsic_mag - 15 + 5*log10(range_) - 2.5*log10(sin(beta) + (pi-beta)*cos(beta))
        return mag

    def sat_location_sun_angle(self, jd: float) -> float:
        """
        Get phase angle between location -- satellite -- sun
        Return value in radians
        """
        sun_rho = sun_pos(jd) - self.location.recef
        sat_rho = self.rho_jd(jd)
        phase_angle = acos(np.dot(sat_rho, sun_rho) / (norm(sat_rho) * norm(sun_rho)))
        return phase_angle

    def determine_visibility(self, d: datetime.datetime) -> Visibility:
        """
        Determine if satellite is visible for single datetime instance
        """
        jd = sum(julian_date_from_datetime(d))
        return self.determine_visibility_jd(jd)

    def determine_visibility_jd(self, jd: float) -> Visibility:
        """
        Determine if satellite is visible for single julian date
        """
        sat_el = self._elevation_at_jd(jd)
        if sat_el < self.aos_at:
            return None
        sun_el = self.location.sun_elevation_jd(jd)
        if sun_el > self.sunrise_dg:
            return Visibility.daylight
        dist = self.satellite.illumination_distance_jd(jd)
        if dist > R_EARTH:
            return Visibility.visible
        return Visibility.unlit
