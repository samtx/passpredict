from __future__ import annotations
from collections import defaultdict
import datetime
import typing
from math import radians, degrees, pi, log10, sin, cos, acos, floor
from functools import cached_property, lru_cache
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.linalg import norm
from orbit_predictor.predictors.pass_iterators import LocationPredictor
from scipy.interpolate import CubicSpline

from . import _rotations
from .time import julian_date_from_datetime
from ._time import jday2datetime
from .solar import sun_pos
from .constants import R_EARTH
from .exceptions import NotReachable, PropagationError
from .locations import Location
from .utils import get_pass_detail_datetime_metadata

if typing.TYPE_CHECKING:
    from .satellites import SatellitePredictor, LLH


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
    ):
        self.aos_dt = aos_dt if aos_dt is not None else None
        self.tca_dt = tca_dt
        self.los_dt = los_dt if los_dt is not None else None
        self.max_elevation = max_elevation
        self.type = type_
        self.vis_begin_dt = vis_begin_dt
        self.vis_end_dt = vis_end_dt

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
    brightness: float = None

    @cached_property
    def midpoint(self):
        """ Return datetime in UTC of midpoint of pass """
        midpt = self.aos.dt + (self.los.dt - self.aos.dt) / 2
        return midpt

    @cached_property
    def duration(self):
        """ Return pass duration in seconds """
        return (self.los.dt - self.aos.dt).total_seconds()

    def __repr__(self):
        return f"<PredictedPass {self.satid} over {repr(self.location)} on {self.aos.dt}"


class Observer(LocationPredictor):
    """
    Predicts passes of a satellite over a given location.
    Exposes an iterable interface.
    Notice that this algorithm is not fully exhaustive,
    see https://github.com/satellogic/orbit-predictor/issues/99 for details.
    """

    def __init__(
        self,
        location: Location,
        satellite: SatellitePredictor,
        max_elevation_gt=0,
        aos_at_dg=0,
        tolerance_s=1.0
    ):
        """
        Initialize Observer but also compute radians for geodetic coordinates
        """
        self.location = location
        self.satellite = satellite
        self.max_elevation_gt = radians(max([max_elevation_gt, aos_at_dg]))
        self.set_minimum_elevation(aos_at_dg)
        self.set_tolerance(tolerance_s)

    @property
    def predictor(self):
        """ For backwards compatibility """
        return self.satellite

    def iter_passes(self, start_date, limit_date=None):
        """Returns one pass each time"""
        current_date = start_date
        while True:
            if self._is_ascending(current_date):
                # we need a descending point
                ascending_date = current_date
                descending_date = self._find_nearest_descending(ascending_date)
                pass_ = self._refine_pass(ascending_date, descending_date)
                if pass_.valid:
                    if limit_date is not None and pass_.aos > limit_date:
                        break
                    predicted_pass = self._build_predicted_pass(pass_)
                    yield predicted_pass
                if limit_date is not None and current_date > limit_date:
                    break
                current_date = pass_.tca + self._orbit_step(0.6)
            else:
                current_date = self._find_nearest_ascending(current_date)

    @property
    def passes_over(self, *a, **kw):
        return self.iter_passes(*a, **kw)

    def get_next_pass(self,
        aos_dt: datetime.datetime,
        *,
        limit_date: datetime.datetime = None,
    ) -> PredictedPass:
        """
        Gets first overpass starting at aos_dt
        """
        pass_ = next(self.iter_passes(aos_dt, limit_date=limit_date))
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

    def set_minimum_elevation(self, elevation: float):
        """  Set minimum elevation for an overpass  """
        self.aos_at = radians(elevation)
        self.aos_at_deg = elevation

    def set_tolerance(self, seconds: float):
        """  Set tolerance variables """
        if seconds <= 0:
            raise Exception("Tolerance must be > 0")
        self.tolerance_s = seconds
        self.tolerance = datetime.timedelta(seconds=seconds)

    def _build_predicted_pass(self, basic_pass: BasicPassInfo):
        """Returns a classic predicted pass"""
        data = defaultdict()
        data.update({
            'satid': self.predictor.satid,
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

    def _find_nearest_descending(self, ascending_date):
        for candidate in self._sample_points(ascending_date):
            if not self._is_ascending(candidate):
                return candidate
        else:
            # logger.error('Could not find a descending pass over %s start date: %s - TLE: %s',
            #              self.location, ascending_date, self.predictor.tle)
            raise PropagationError("Can not find an descending phase")

    def _find_nearest_ascending(self, descending_date):
        for candidate in self._sample_points(descending_date):
            if self._is_ascending(candidate):
                return candidate
        else:
            # logger.error('Could not find an ascending pass over %s start date: %s - TLE: %s',
            #              self.location, descending_date, self.predictor.tle)
            raise PropagationError('Can not find an ascending phase')

    def _sample_points(self, date):
        """Helper method to found ascending or descending phases of elevation"""
        start = date
        end = date + self._orbit_step(0.99)
        mid = self.midpoint(start, end)
        mid_right = self.midpoint(mid, end)
        mid_left = self.midpoint(start, mid)
        return [end, mid, mid_right, mid_left]

    def _refine_pass(self, ascending_date, descending_date) -> BasicPassInfo:
        tca_dt = self._find_tca(ascending_date, descending_date)
        elevation = self._elevation_at(tca_dt)
        if elevation > self.max_elevation_gt:
            aos_dt = self._find_aos(tca_dt)
            los_dt = self._find_los(tca_dt)
        else:
            aos_dt = los_dt = None
        return BasicPassInfo(aos_dt, tca_dt, los_dt, elevation)
        # Find visual pass details
        # First, get endpoints of when location is not sunlit
        # Use cubic splines to find sun elevation
        jd0 = sum(julian_date_from_datetime(aos_dt))
        jdf = sum(julian_date_from_datetime(los_dt))
        jd = np.linspace(jd0, jdf, 5)  # use 5 points for spline
        sun_el_fn = lambda j: self.location.sun_elevation_jd(j) + 6
        el = np.array([sun_el_fn(j) for j in jd])
        if np.min(el) > 0:
            # entire pass in sunlit
            return BasicPassInfo(aos_dt, tca_dt, los_dt, elevation, type_=PassType.daylight)

        tol = 1/86400  # one second
        # part of the pass is in darkness. Find new jd0, jdf
        sun_el = CubicSpline(jd, el, bc_type='natural')
        for root in sun_el.roots(extrapolate=False):
            tmp1 = sun_el(root - tol)
            tmp2 = sun_el(root + tol)
            if tmp1 < tmp2:
                # sun elevation is decreasing
                jd0 = root
            else:
                jdf = root
        # Now use jd0 and jdf to find when satellite is illuminated by sun
        jd = np.linspace(jd0, jdf, 5)  # use 10 points for spline
        illum_fn = lambda j: self.satellite.illumination_distance_jd(j) - R_EARTH
        illum_pts = np.array([illum_fn(j) for j in jd])
        if np.max(illum_pts) < 0:
            # entire pass is in shadow
            return BasicPassInfo(aos_dt, tca_dt, los_dt, elevation, type_=PassType.unlit)

        # part of the pass is visible
        illum = CubicSpline(jd, illum_pts, bc_type='natural')
        for root in illum.roots(extrapolate=False):
            tmp1 = illum(root - tol)
            tmp2 = illum(root + tol)
            if tmp1 < tmp2:
                # satellite is going into shadow
                jdf = root
            else:
                # satellite is coming out of shadow
                jd0 = root
        # Set visible start and end points for Pass
        vis_begin_dt = jday2datetime(jd0)
        vis_end_dt = jday2datetime(jdf)
        return BasicPassInfo(
            aos_dt, tca_dt, los_dt, elevation, type_=PassType.visible,
            vis_begin_dt=vis_begin_dt, vis_end_dt=vis_end_dt,
        )

    def _find_tca(self, ascending_date, descending_date):
        while not self._precision_reached(ascending_date, descending_date):
            midpoint = self.midpoint(ascending_date, descending_date)
            if self._is_ascending(midpoint):
                ascending_date = midpoint
            else:
                descending_date = midpoint
        return ascending_date

    def _precision_reached(self, start, end):
        return end - start <= self.tolerance

    @staticmethod
    def midpoint(start, end):
        """Returns the midpoint between two dates"""
        return start + (end - start) / 2

    @lru_cache(maxsize=128)
    def _elevation_at(self, when_utc):
        position = self.predictor.get_only_position(when_utc)
        return self.location.elevation_for(position)

    def _is_ascending(self, when_utc):
        """Check is elevation is ascending or descending on a given point"""
        elevation = self._elevation_at(when_utc)
        next_elevation = self._elevation_at(when_utc + self.tolerance)
        return elevation <= next_elevation

    def _orbit_step(self, size):
        """Returns a time step, that will make the satellite advance a given number of orbits"""
        step_in_radians = size * 2 * pi
        seconds = (step_in_radians / self.predictor.mean_motion) * 60
        return datetime.timedelta(seconds=seconds)

    def _find_aos(self, tca):
        end = tca
        start = tca - self._orbit_step(0.34)  # On third of the orbit
        elevation = self._elevation_at(start)
        assert elevation < 0
        while not self._precision_reached(start, end):
            midpoint = self.midpoint(start, end)
            elevation = self._elevation_at(midpoint)
            if elevation < self.aos_at:
                start = midpoint
            else:
                end = midpoint
        return end

    def _find_los(self, tca):
        start = tca
        end = tca + self._orbit_step(0.34)
        while not self._precision_reached(start, end):
            midpoint = self.midpoint(start, end)
            elevation = self._elevation_at(midpoint)
            if elevation < self.aos_at:
                end = midpoint
            else:
                start = midpoint
        return start

    def razel(self, datetime: datetime.datetime) -> RangeAzEl:
        """
        Get range, azimuth, and elevation for datetime
        """
        satellite_ecef = self.predictor.get_only_position(datetime)
        range_, az, el = _rotations.razel(
            self.location.latitude_rad, self.location.longitude_rad, self.location.recef, satellite_ecef
        )
        return RangeAzEl(range_, az, el)

    @lru_cache(maxsize=128)
    def range_at_jd(self, jd: float) -> float:
        """
        Get slant range magnitude from location to satellite [km]
        """
        satellite_ecef = self.predictor.get_only_position_jd(jd)
        range_ = _rotations.range_at(
            self.location.latitude_rad,
            self.location.longitude_rad,
            self.location.recef,
            satellite_ecef,
        )
        return range_

    def point(self, datetime: datetime.datetime) -> PassPoint:
        """
        Get PassPoint with range, azimuth, and elevation data for datetime
        """
        rnazel = self.razel(datetime)
        pt = PassPoint(datetime, rnazel.range, rnazel.az, rnazel.el)
        return pt

    @lru_cache(maxsize=128)
    def rho_jd(self, jd: float) -> np.ndarray:
        """
        Get topocentric ECEF vector from location to satellite
        """
        rho = np.empty(3, dtype=np.double)
        sat_recef = self.predictor.get_only_position_jd(jd)
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
