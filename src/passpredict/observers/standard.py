from __future__ import annotations
from collections import defaultdict
import datetime
import typing
from math import radians, pi

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import root_scalar

from .base import ObserverBase, PredictedPass, BasicPassInfo, PassType, Visibility
from .functions import find_root, julian_date_sum
from ..time import julian_date_from_datetime
from .._time import jday2datetime_us
from ..constants import R_EARTH
from ..exceptions import PropagationError

if typing.TYPE_CHECKING:
    from ..satellites import SatellitePredictor, LLH
    from ..locations import Location


def _make_utc(d: datetime.datetime) -> datetime.datetime:
    """ Make a datetime a UTC timezone aware datetime """
    if not d:
        return d
    if d.tzinfo:
        d = d.astimezone(datetime.timezone.utc)
    else:
        d = d.replace(tzinfo=datetime.timezone.utc)
    return d


class Observer(ObserverBase):
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
        **kwargs,
    ):
        """
        Initialize Observer
        """
        super().__init__(location, satellite, **kwargs)

    def iter_passes(self, start_date, limit_date=None, visible_only=False):
        """Returns one pass each time"""
        start_date = _make_utc(start_date)
        limit_date = _make_utc(limit_date)
        current_date = start_date
        while True:
            if self._is_ascending(current_date):
                # we need a descending point
                ascending_date = current_date
                descending_date = self._find_nearest_descending(ascending_date)
                pass_ = self._refine_pass(ascending_date, descending_date)
                if self._is_pass_valid(pass_, visible_only=visible_only):
                    if limit_date is not None and pass_.aos > limit_date:
                        break
                    predicted_pass = self._build_predicted_pass(pass_)
                    yield predicted_pass
                if limit_date is not None and current_date > limit_date:
                    break
                current_date = pass_.tca + self._orbit_step(0.6)
            else:
                current_date = self._find_nearest_ascending(current_date)

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
        tz = basic_pass.aos.tzinfo
        if basic_pass.vis_begin_dt:
            data['vis_begin'] = self.point(basic_pass.vis_begin_dt.astimezone(tz))
        if basic_pass.vis_end_dt:
            data['vis_end'] = self.point(basic_pass.vis_end_dt.astimezone(tz))
        if basic_pass.vis_tca_dt:
            data['vis_tca'] = self.point(basic_pass.vis_tca_dt.astimezone(tz))
        if (basic_pass.type == Visibility.visible) and (self.satellite.intrinsic_mag is not None):
            brightness = 999
            for pt in ('vis_begin', 'vis_tca', 'vis_end',):
                b = data[pt].brightness
                if b and b < brightness:
                    brightness = b
            data['brightness'] = brightness
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
        jd0 = sum(julian_date_from_datetime(aos_dt))
        jdf = sum(julian_date_from_datetime(los_dt))
        jd = np.linspace(jd0, jdf, 5)  # use 5 points for spline
        sun_el_fn = lambda j: self.location.sun_elevation_jd(j) - self.sunrise_dg
        el = np.array([sun_el_fn(j) for j in jd])
        if np.min(el) > 0:
            # entire pass in sunlit
            return BasicPassInfo(aos_dt, tca_dt, los_dt, elevation, type_=PassType.daylight)
        # part of the pass is in darkness.
        if el[0]*el[-1] < 0:
            # only part of the pass is sunlit. Find new jd0, jdf
            result = root_scalar(sun_el_fn, method='bisect', bracket=(jd0, jdf), x0=jd0, xtol=self.jd_tol)
            x = result.root
            tmp1 = sun_el_fn(x - self.jd_tol)
            tmp2 = sun_el_fn(x + self.jd_tol)
            if tmp1 < tmp2:
                # sun elevation is decreasing
                jd0 = x
            else:
                jdf = x

        # Now use jd0 and jdf to find when satellite is illuminated by sun
        jd = np.linspace(jd0, jdf, 5)  # use 5 points for spline
        illum_fn = lambda j: self.satellite.illumination_distance_jd(j) - R_EARTH
        illum_pts = np.array([illum_fn(j) for j in jd])
        if np.max(illum_pts) < 0:
            # entire pass is in shadow
            return BasicPassInfo(aos_dt, tca_dt, los_dt, elevation, type_=PassType.unlit)

        if illum_pts[0]*illum_pts[-1] < 0:
            # the satellite is visible for only part of the pass. Find new jd0, jdf
            result = root_scalar(illum_fn, method='bisect', bracket=(jd0, jdf), x0=jd0, xtol=self.jd_tol)
            x = result.root
            tmp1 = illum_fn(x - self.jd_tol)
            tmp2 = illum_fn(x + self.jd_tol)
            if tmp1 < tmp2:
                # satellite is coming out of shadow
                jd0 = x
            else:
                # satellite is going into shadow
                jdf = x
        # Set visible start and end points for Pass
        vis_begin_dt = jday2datetime_us(jd0)
        vis_end_dt = jday2datetime_us(jdf)
        # Find maximum elevation during visible period
        if vis_begin_dt <= tca_dt <= vis_end_dt:
            vis_tca_dt = tca_dt
        elif self._elevation_at(vis_begin_dt) > self._elevation_at(vis_end_dt):
            vis_tca_dt = vis_begin_dt
        else:
            vis_tca_dt = vis_end_dt
        return BasicPassInfo(
            aos_dt, tca_dt, los_dt, elevation, type_=PassType.visible,
            vis_begin_dt=vis_begin_dt, vis_end_dt=vis_end_dt, vis_tca_dt=vis_tca_dt,
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
        # assert elevation < 0
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
