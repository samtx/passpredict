from __future__ import annotations
import datetime as dt
from typing import NamedTuple
from math import radians, degrees, pi
from functools import cached_property
from dataclasses import dataclass
import logging
from enum import Enum

import numpy as np
from orbit_predictor.predictors.pass_iterators import LocationPredictor as BaseLocationPredictor
from orbit_predictor.predictors.pass_iterators import PredictedPass as BasePredictedPass
from orbit_predictor.predictors.accurate import HighAccuracyTLEPredictor
from orbit_predictor.locations import Location

from . import _rotations
from .exceptions import NotReachable, PropagationError
from .locations import Location
from .utils import direction_from_azimuth
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .sources import TLESource, TLE


class RangeAzEl(NamedTuple):
    range: float  # km
    az: float     # deg
    el: float     # deg


class PassType(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'


@dataclass
class PassPoint:
    dt: dt.datetime    # datetime UTC
    range_: float      # km
    azimuth: float     # deg
    elevation: float   # deg

    @property
    def range(self):
        return self.range_

    @cached_property
    def direction(self) -> str:
        ''' Return direction from azimuth degree '''
        return direction_from_azimuth(self.azimuth)


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
        aos_dt: dt.datetime,
        tca_dt: dt.datetime,
        los_dt: dt.datetime,
        max_elevation: float
    ):
        self.aos_dt = aos_dt if aos_dt is not None else None
        self.tca_dt = tca_dt
        self.los_dt = los_dt if los_dt is not None else None
        self.max_elevation = max_elevation

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
    def duration(self) -> dt.timedelta:
        return self.los - self.aos


@dataclass(frozen=True)
class PredictedPass:
    satid: int
    location: Location
    aos: PassPoint
    tca: PassPoint
    los: PassPoint
    brightness: float = None
    type: PassType = None

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


class SatellitePredictor(HighAccuracyTLEPredictor):
    """
    Predictor for satellite overpasses. Uses sgp4 for propagation
    """
    def __init__(self, satid: int, source: TLESource = None):
        """
        Params:
            satid: int = NORAD id for satellite
            source: PasspredictTLESource
        """
        self.satid = satid
        self._source = source
        if self._source:
            self.tle = self._source.get_tle(satid, dt.datetime.utcnow())
            self._propagator = self._get_propagator()
        else:
            self.tle = None
            self._propagator = None

    def get_tle(self):
        self.tle = self._source.get_tle(self.sate_id, dt.datetime.now(dt.timezone.utc))

    @classmethod
    def from_tle(cls, tle: TLE) -> SatellitePredictor:
        predictor = cls(tle.satid)
        predictor.tle = tle
        predictor.set_propagator()
        return predictor

    @property
    def sate_id(self):
        return self.satid

    def set_propagator(self):
        self._propagator = self._get_propagator()

    def get_only_position(self, datetime: dt.datetime) -> np.ndarray:
        """
        Get satellite position in ECEF coordinates [km]
        """
        pos_tuple = super().get_only_position(datetime)
        return np.array(pos_tuple)

    @property
    def passes_over(self, *a, **kw):
        return self.pass_iterator(*a, **kw)

    def pass_iterator(
        self,
        location: Location,
        when_utc: dt.datetime,
        *,
        limit_date: dt.datetime = None,
        max_elevation_gt: float = 0,
        aos_at_dg: float = 0,
        tolerance_s: float = 1.0,
    ):
        return LocationPredictor(
            location,
            self,
            when_utc,
            limit_date,
            max_elevation_gt,
            aos_at_dg,
            tolerance_s=tolerance_s
        )

    def get_next_pass(self,
        location: Location,
        aos_dt: dt.datetime = None,
        *,
        max_elevation_gt: float = 5,
        aos_at_dg: float = 0,
        limit_date: dt.datetime = None,
        tolerance_s: float = 1
        ) -> PredictedPass:
        """
        Gets first overpass starting at aos_dt
        """
        if aos_dt is None:
            aos_dt = dt.datetime.utcnow()
        for pass_ in self.pass_iterator(location, aos_dt, limit_date=limit_date, max_elevation_gt=max_elevation_gt, aos_at_dg=aos_at_dg, tolerance_s=tolerance_s):
            return pass_
        else:
            raise NotReachable('Propagation limit date exceeded')


class LocationPredictor(BaseLocationPredictor):
    """Predicts passes over a given location.
    Exposes an iterable interface.
    Notice that this algorithm is not fully exhaustive,
    see https://github.com/satellogic/orbit-predictor/issues/99 for details.
    """

    def __init__(self, *a, **kw):
        """
        Initialize LocationPredictor but also compute radians for geodetic coordinates
        """
        super().__init__(*a, **kw)
        self.aos_at_deg = degrees(self.aos_at)  # use degrees instead of radians
        self.location_lat_rad = radians(self.location.latitude_deg)
        self.location_lon_rad = radians(self.location.longitude_deg)
        self.location_ecef = np.array(self.location.position_ecef)

    def iter_passes(self):
        """Returns one pass each time"""
        current_date = self.start_date
        while True:
            if self._is_ascending(current_date):
                # we need a descending point
                ascending_date = current_date
                descending_date = self._find_nearest_descending(ascending_date)
                pass_ = self._refine_pass(ascending_date, descending_date)
                if pass_.valid:
                    if self.limit_date is not None and pass_.aos > self.limit_date:
                        break
                    # if pass_.max_elevation_deg < self.max_elevation_gt:
                    #     break
                    predicted_pass = self._build_predicted_pass(pass_)
                    # if (predicted_pass.aos.elevation < self.aos_at_deg) or (predicted_pass.los.elevation < self.aos_at_deg):
                    #     break
                    yield predicted_pass

                if self.limit_date is not None and current_date > self.limit_date:
                    break

                current_date = pass_.tca + self._orbit_step(0.6)

            else:
                current_date = self._find_nearest_ascending(current_date)

    def _build_predicted_pass(self, basic_pass: BasicPassInfo):
        """Returns a classic predicted pass"""
        aos = self.point(basic_pass.aos_dt)
        tca = self.point(basic_pass.tca_dt)
        los = self.point(basic_pass.los_dt)
        return PredictedPass(self.predictor.satid, self.location, aos, tca, los)

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

    def _refine_pass(self, ascending_date, descending_date):
        tca_dt = self._find_tca(ascending_date, descending_date)
        elevation = self._elevation_at(tca_dt)
        if elevation > self.max_elevation_gt:
            aos_dt = self._find_aos(tca_dt)
            los_dt = self._find_los(tca_dt)
        else:
            aos_dt = los_dt = None
        return BasicPassInfo(aos_dt, tca_dt, los_dt, elevation)

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
        return dt.timedelta(seconds=seconds)

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

    def razel(self, datetime: dt.datetime) -> RangeAzEl:
        """
        Get range, azimuth, and elevation for datetime
        """
        satellite_ecef = self.predictor.get_only_position(datetime)
        range_, az, el = _rotations.razel(
            self.location_lat_rad, self.location_lon_rad, self.location_ecef, satellite_ecef
        )
        return RangeAzEl(range_, az, el)

    def point(self, datetime: dt.datetime) -> PassPoint:
        """
        Get PassPoint with range, azimuth, and elevation data for datetime
        """
        rnazel = self.razel(datetime)
        pt = PassPoint(datetime, rnazel.range, rnazel.az, rnazel.el)
        return pt