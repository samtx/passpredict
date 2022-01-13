from __future__ import annotations
import datetime
from math import pi, log10, sin, cos, acos
import typing
from functools import lru_cache
from collections import defaultdict
from abc import abstractmethod

import numpy as np
from numpy.linalg import norm
from orbit_predictor.predictors.pass_iterators import LocationPredictor

from .predicted_pass import RangeAzEl, BasicPassInfo, PassPoint, PredictedPass
from .. import _rotations
from ..utils import get_pass_detail_datetime_metadata
from ..exceptions import NotReachable
from ..solar import sun_pos

if typing.TYPE_CHECKING:
    from ..satellites import LLH


class ObserverBase(LocationPredictor):

    @property
    def predictor(self):
        """ For backwards compatibility with orbit-predictor location predictor """
        return self.satellite

    @abstractmethod
    def iter_passes(self, start_date, limit_date=None):
        raise NotImplementedError

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

    def point(self, datetime: datetime.datetime) -> PassPoint:
        """
        Get PassPoint with range, azimuth, and elevation data for datetime
        """
        rnazel = self.razel(datetime)
        pt = PassPoint(datetime, rnazel.range, rnazel.az, rnazel.el)
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