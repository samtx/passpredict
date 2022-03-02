from __future__ import annotations
import datetime
from typing import TYPE_CHECKING, List, Tuple
import warnings
from functools import lru_cache
from collections import defaultdict
from math import sin, cos, log10, acos, pi, radians, degrees

import numpy as np
from numpy.linalg import norm

from .core import (
    PredictedPass,
    PassType,
    Visibility,
    BasicPassInfo,
    PassPoint,
    RangeAzEl,
)
from .functions import make_utc
from .orbit_predictor import orbit_predictor_iterator
from .brute_force import brute_force_iterator
from .. import _rotations
from .._time import mjd2datetime, datetime2mjd
from ..solar import sun_pos_mjd
from ..constants import R_EARTH
from ..utils import get_pass_detail_datetime_metadata

if TYPE_CHECKING:
    from ..locations import Location
    from ..satellites import LLH


class Observer:

    def __init__(
        self,
        location: Location,
        satellite,
    ):
        self.location = location
        self.satellite = satellite

    @property
    def predictor(self):
        """ For backwards compatibility with orbit-predictor location predictor """
        warnings.warn(".predictor is deprecated; use .satellite", DeprecationWarning)
        return self.satellite

    def iter_passes(
        self,
        start_date: datetime.datetime,
        limit_date: datetime.datetime = None,
        method: str = 'op',
        visible_only: bool = False,
        **options
    ):
        start_mjd = datetime2mjd(make_utc(start_date))

        # Get options with defaults
        aos_at_dg = options.get('aos_at_dg', 0)
        max_elevation_gt = options.get('max_elevation_gt', 0)

        if 'tolerance_s' in options:
            warnings.warn("tolerance_s has been replaced with tol", DeprecationWarning)
            if 'tol' not in options:
                options['tol'] = options['tolerance_s']
        tol = options.get('tol', 1)   # in seconds
        if tol < 0:
            raise Exception("tol must be > 0")
        sunrise_dg = options.get('sunrise_dg', -6)

        # Derived options
        aos_at = radians(aos_at_dg)

        if method == 'op':
            for pass_ in orbit_predictor_iterator(
                self, start_date, limit_date, aos_at_dg=aos_at_dg,
                max_elevation_gt=max_elevation_gt, tol=tol,
                sunrise_dg=sunrise_dg,
            ):
                if _is_valid(pass_, visible_only, max_elevation_gt, start_mjd):
                    predicted_pass = self._build_predicted_pass(
                        pass_, aos_at=aos_at, sunrise_dg=sunrise_dg
                    )
                    yield predicted_pass

        elif method == 'brute':
            # Get options for brute force iterator
            time_step = options.get('time_step', 10)  # in seconds
            if time_step <= 0:
                raise Exception("time_step must be > 0")
            for pass_ in brute_force_iterator(
                self, start_date, limit_date, aos_at_dg=aos_at_dg,
                max_elevation_gt=max_elevation_gt, tol=tol,
                sunrise_dg=sunrise_dg, time_step=time_step,
            ):
                if _is_valid(pass_, visible_only, max_elevation_gt, start_mjd):
                    predicted_pass = self._build_predicted_pass(
                        pass_, aos_at=aos_at, sunrise_dg=sunrise_dg
                    )
                    yield predicted_pass

    def pass_list(
        self,
        start_date: datetime.datetime,
        limit_date: datetime.datetime,
        *,
        visible_only: bool = False,
        **kwargs,
    ) -> List[PredictedPass]:
        """
        Compute overpass predictions. Return list of PredictedPass objects

        Parameters
        ----------
        start_date : datetime.datetime
            The start date of the overpass predictions
        limit_date : datetime.datetime
            The end date of the overpass predictions
        visible_only : bool, optional
            Return only visible overpasses, by default False

        Returns
        -------
        List[PredictedPass]
            A list of PredictedPass objects
        """
        return list(self.iter_passes(
            start_date, limit_date=limit_date,
            visible_only=visible_only, **kwargs,
        ))

    def get_next_pass(
        self,
        *a,
        **kw,
    ) -> PredictedPass:
        """ Gets first overpass starting at aos_dt """
        warnings.warn("get_next_pass() is deprecated; use next_pass().", DeprecationWarning)
        return self.next_pass(*a, **kw)

    def next_pass(
        self,
        aos_dt: datetime.datetime,
        *,
        limit_date: datetime.datetime = None,
        visible_only: bool = False,
        **kw,
    ) -> PredictedPass:
        """
        Gets first overpass starting at aos_dt
        """
        pass_ = next(
            self.iter_passes(aos_dt, limit_date=limit_date, visible_only=visible_only, **kw),
            None
        )
        return pass_

    def next_pass_detail(
        self,
        aos_dt: datetime.datetime,
        *,
        limit_date: datetime.datetime = None,
        delta_s: float = 10,
        pad_minutes: int = 5,
    ) -> Tuple[PredictedPass, LLH]:
        """
        Add details to PredictedPass
        Evaluate position and velocity properties for each delta_s seconds
        """
        pass_ = self.next_pass(aos_dt, limit_date=limit_date)
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

    def _build_predicted_pass(
        self,
        basic_pass: BasicPassInfo,
        aos_at: float,
        sunrise_dg: float
    ):
        """Returns a classic predicted pass"""
        data = defaultdict()
        kw = {
            'aos_at': aos_at,
            'sunrise_dg': sunrise_dg,
        }
        data.update({
            'satid': self.satellite.satid,
            'location': self.location,
            'type': basic_pass.type,
            'aos': self._point_mjd(basic_pass.aos_mjd, **kw),
            'tca': self._point_mjd(basic_pass.tca_mjd, **kw),
            'los': self._point_mjd(basic_pass.los_mjd, **kw),
        })
        if basic_pass.vis_begin_mjd:
            data['vis_begin'] = self._point_mjd(basic_pass.vis_begin_mjd, **kw)
        if basic_pass.vis_end_mjd:
            data['vis_end'] = self._point_mjd(basic_pass.vis_end_mjd, **kw)
        if basic_pass.vis_tca_mjd:
            data['vis_tca'] = self._point_mjd(basic_pass.vis_tca_mjd, **kw)
        if (basic_pass.type == Visibility.visible) and (self.satellite.intrinsic_mag is not None):
            brightness = 999
            for pt in ('vis_begin', 'vis_tca', 'vis_end',):
                b = data[pt].brightness
                if b and b < brightness:
                    brightness = b
            data['brightness'] = brightness
        return PredictedPass(**data)

    def elevation(self, d: datetime.datetime):
        """  Return elevation of object in degrees  """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        el = self._elevation_mjd(mjd)
        return degrees(el)

    @lru_cache(maxsize=16)
    def _elevation_mjd(self, mjd: float) -> float:
        """  Returns elevation of object in radians  """
        sat_recef = self.satellite._position_ecef_mjd(mjd)
        coslatcoslon, coslatsinlon, sinlat = self.location._cached_elevation_calculation_data
        return _rotations.elevation_at_rad(coslatcoslon, coslatsinlon, sinlat, self.location.recef, sat_recef)

    def razel(self, d: datetime.datetime) -> RangeAzEl:
        """
        Get range, azimuth, and elevation for datetime
        """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        return self._razel_mjd(mjd)

    def _razel_mjd(self, mjd: float) -> RangeAzEl:
        """
        Get range, azimuth, and elevation for mjd time
        """
        satellite_ecef = self.satellite._position_ecef_mjd(mjd)
        range_, az, el = _rotations.razel(
            self.location.latitude_rad, self.location.longitude_rad, self.location.recef, satellite_ecef
        )
        return RangeAzEl(range_, az, el)

    def point(self, d: datetime.datetime, **kw) -> PassPoint:
        """
        Get PassPoint with range, azimuth, and elevation data for mjd time
        """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        return self._point_mjd(mjd)

    def _point_mjd(self, mjd: float, **kw) -> PassPoint:
        """
        Get PassPoint with range, azimuth, and elevation data for mjd time
        """
        rnazel = self._razel_mjd(mjd)
        vis_state = self._determine_visibility_mjd(mjd, **kw)
        if (vis_state == Visibility.visible) and (self.satellite.intrinsic_mag is not None):
            brightness = self._brightness_mjd(mjd)
        else:
            brightness = None
        d = mjd2datetime(mjd)
        pt = PassPoint(
            d, rnazel.range, rnazel.az, rnazel.el, type=vis_state, brightness=brightness
        )
        return pt

    def rho(self, d: datetime.datetime) -> np.ndarray:
        """
        Get topocentric ECEF vector from location to satellite
        """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        return self._rho_mjd(mjd)

    @lru_cache(maxsize=16)
    def _rho_mjd(self, mjd: float) -> np.ndarray:
        """
        Get topocentric ECEF vector from location to satellite
        """
        rho = np.empty(3, dtype=np.double)
        sat_recef = self.satellite._position_ecef_mjd(mjd)
        _rotations.ecef_to_rhosez(self.location.latitude_rad, self.location.longitude_rad, self.location.recef, sat_recef, rho)
        return rho

    def _brightness_mjd(self, mjd: float) -> float:
        beta, range_ = self._sat_location_sun_angle_mjd(mjd)
        beta = pi - beta
        mag = self.satellite.intrinsic_mag - 15 + 5*log10(range_) - 2.5*log10(sin(beta) + (pi-beta)*cos(beta))
        return mag

    def brightness(self, d: datetime.datetime) -> float:
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        return self._brightness_mjd(mjd)

    def _sat_location_sun_angle_mjd(self, mjd: float):
        """
        Get phase angle between location -- satellite -- sun
        Return value in radians
        Also return slant range to satellite
        """
        sun_rho = sun_pos_mjd(mjd) - self.location.recef
        sat_rho = self._rho_mjd(mjd)
        range_ = norm(sat_rho)
        phase_angle = acos(np.dot(sat_rho, sun_rho) / (range_ * norm(sun_rho)))
        return phase_angle, range_

    def sat_location_sun_angle(self, d: datetime.datetime) -> float:
        """
        Get phase angle between location -- satellite -- sun
        Return value in radians
        """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        phase_angle, _ = self._sat_location_sun_angle_mjd(mjd)
        return phase_angle

    def determine_visibility(self, d: datetime.datetime, **kw) -> Visibility:
        """
        Determine if satellite is visible for single datetime instance
        """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        return self._determine_visibility_mjd(mjd, **kw)

    def _determine_visibility_mjd(
        self,
        mjd: float,
        *,
        aos_at: float = 0,
        sunrise_dg: float = -6,
        **kw
    ) -> Visibility:
        """
        Determine if satellite is visible for single modified julian date
        """
        sat_el = self._elevation_mjd(mjd)
        if sat_el < aos_at:
            return None
        sun_el = self.location._sun_elevation_mjd(mjd)
        if sun_el > sunrise_dg:
            return Visibility.daylight
        dist = self.satellite._illumination_distance_mjd(mjd)
        if dist > R_EARTH:
            return Visibility.visible
        return Visibility.unlit


def _is_valid(pass_, visible_only, max_elevation_gt, start_mjd):
    if not pass_:
        return False
    if pass_.aos_mjd < start_mjd:
        return False
    if visible_only and pass_.type != PassType.visible:
        return False
    return pass_.max_elevation >= max_elevation_gt
