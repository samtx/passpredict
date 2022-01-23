from __future__ import annotations
import datetime
import typing
from functools import lru_cache

import numpy as np
from orbit_predictor.predictors.accurate import HighAccuracyTLEPredictor
from orbit_predictor.coordinate_systems import ecef_to_llh, eci_to_ecef
from sgp4.api import SGP4_ERRORS, Satrec
from sgp4.model import WGS84
from sgp4.propagation import gstime

from ..time import julian_date_from_datetime
from ..solar import sun_pos
from .. import _solar
from ..exceptions import PropagationError
from ..constants import R_EARTH

if typing.TYPE_CHECKING:
    from ..sources import PasspredictTLESource, TLE


class LLH(typing.NamedTuple):
    datetime: typing.Sequence[datetime.datetime]
    latitude: np.typing.NDArray[float]
    longitude: np.typing.NDArray[float]
    altitude: np.typing.NDArray[float]


class SatellitePredictorBase(HighAccuracyTLEPredictor):
    """
    Predictor for satellite overpasses.
    """
    def __init__(self, satid: int, source: PasspredictTLESource = None):
        """
        Params:
            satid: int = NORAD id for satellite
            source: PasspredictTLESource
        """
        self.satid = satid
        if source:
            self.tle = source.get_tle(satid)
            self.name = self.tle.name
            self._propagator = self.get_propagator(self.tle.lines)
        self.intrinsic_mag = 1.0   # ISS is -1.8

    @classmethod
    def from_tle(cls, tle: TLE) -> SatellitePredictorBase:
        sat = cls(tle.satid)
        sat.tle = tle
        sat.name = tle.name
        sat._propagator = sat.get_propagator(tle.lines)
        return sat

    @property
    def sate_id(self):
        return self.satid

    def __repr__(self):
        return f"<{self.__class__.__name__} satid={self.satid} (TLE epoch {self.tle.epoch})>"

    def get_propagator(self, lines):
        tle_line_1, tle_line_2 = lines
        return Satrec.twoline2rv(tle_line_1, tle_line_2, WGS84)

    def get_only_position(self, datetime: datetime.datetime) -> np.ndarray:
        """
        Get satellite position in ECEF coordinates [km]
        """
        pos_tuple = super().get_only_position(datetime)
        return np.array(pos_tuple)

    @lru_cache(maxsize=1800)  # Max cache, 30 minutes
    def get_only_position_jd(self, jd: float) -> np.ndarray:
        """
        Get satellite position in ECEF coordinates [km]
        """
        status, position_eci, _ = self._propagator.sgp4(jd, 0.0)
        if status != 0:
            raise PropagationError(SGP4_ERRORS[status])
        gmst = gstime(jd)
        pos_tuple = eci_to_ecef(position_eci, gmst)
        return np.array(pos_tuple)

    def get_llh(self, datetime: datetime.datetime) -> np.ndarray:
        raise NotImplementedError

    def get_position_detail(self, start_date, n_steps, time_step):
        """
        Get satellite subpoints and details
        """
        latitude = np.empty(n_steps)
        longitude = np.empty(n_steps)
        altitude = np.empty(n_steps)
        datetime = [None] * n_steps
        dt = start_date
        for i in range(n_steps):
            recef = self.get_only_position(dt)
            lat, lon, h = ecef_to_llh(recef)
            latitude[i] = lat
            longitude[i] = lon
            altitude[i] = h
            datetime[i] = dt
            dt += time_step
        return LLH(datetime, latitude, longitude, altitude)

    def is_illuminated(self, d: datetime.datetime) -> bool:
        jd = julian_date_from_datetime(d)
        rsun = sun_pos(jd)
        rsat = self.get_only_position(d)
        is_illum = _solar.is_sat_illuminated(rsat, rsun)
        return is_illum

    def illumination_distance_jd(self, jd: float) -> float:
        rsun = sun_pos(jd)
        rsat = self.get_only_position_jd(jd)
        dist = _solar.sat_illumination_distance(rsat, rsun)
        return dist