from __future__ import annotations
import datetime
import typing
import abc

import numpy as np
from orbit_predictor.coordinate_systems import ecef_to_llh

from ..time import make_utc
from .._time import datetime2mjd
from .. import _solar
from ..constants import R_EARTH

if typing.TYPE_CHECKING:
    from ..sources import PasspredictTLESource, TLE


class LLH(typing.NamedTuple):
    datetime: typing.Sequence[datetime.datetime]
    latitude: np.typing.NDArray[float]
    longitude: np.typing.NDArray[float]
    altitude: np.typing.NDArray[float]


class SatellitePropagatorBase(abc.ABC):
    """
    Propagator for satellite position and velocity.
    """

    @property
    def sate_id(self):
        return self.satid

    def __repr__(self):
        return f"<{self.__class__.__name__} satid={self.satid} (TLE epoch {self.tle.epoch})>"

    @abc.abstractproperty
    def mean_motion(self):
        """  Mean motion, in radians per minute  """
        raise NotImplementedError

    def get_only_position(self, d: datetime.datetime) -> np.ndarray:
        """
        Get satellite position in ECEF coordinates [km]
        """
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        return self._position_ecef_mjd(mjd)

    @abc.abstractmethod
    def _position_ecef_mjd(self, mjd: float) -> np.ndarray:
        """ Use modified Julian date """
        raise NotImplementedError

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
        d2 = make_utc(d)
        mjd = datetime2mjd(d2)
        dist = self._illumination_distance_mjd(mjd)
        return dist > R_EARTH

    def _illumination_distance_mjd(self, mjd: float) -> float:
        rsun = _solar.sun_pos_mjd(mjd)
        rsat = self._position_ecef_mjd(mjd)
        dist = _solar.sat_illumination_distance(rsat, rsun)
        return dist

