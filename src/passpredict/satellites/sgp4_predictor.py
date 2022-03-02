from __future__ import annotations
import typing
from functools import cached_property, lru_cache
from math import cos

import numpy as np
from sgp4.api import SGP4_ERRORS, Satrec
from sgp4.model import WGS84
from sgp4.earth_gravity import wgs84 as wgs84_params

from .base import SatellitePropagatorBase
from .._time import mjd2jdfr
from .._rotations import teme2ecef
from ..exceptions import PropagationError

if typing.TYPE_CHECKING:
    from ..sources import TLE


class SGP4Propagator(SatellitePropagatorBase):
    """
    SGP4 Propagator for satellite position and velocity.
    """
    def __init__(self, *a, **kw):
        """
        Params:
            satid: int = NORAD id for satellite
            source: PasspredictTLESource
        """
        super().__init__(*a, **kw)
        if self.tle:
            self.set_propagator()

    @classmethod
    def from_tle(cls, tle: TLE) -> SGP4Propagator:
        sat = cls(tle.satid)
        sat.tle = tle
        sat.name = tle.name
        sat.set_propagator()
        return sat

    def set_propagator(self):
        self._propagator = Satrec.twoline2rv(
            self.tle.lines[0], self.tle.lines[1], WGS84
        )

    @cached_property
    def mean_motion(self):
        """  Mean motion, in radians per minute  """
        return unkozai(
            self._propagator.no_kozai,
            self._propagator.ecco,
            self._propagator.inclo,
            wgs84_params
        )

    @lru_cache(maxsize=1200)
    def _position_ecef_mjd(self, mjd: float) -> np.ndarray:
        """ Use modified Julian date """
        jd, jdfr = mjd2jdfr(mjd)
        status, position_eci, _ = self._propagator.sgp4(jd, jdfr)
        if status != 0:
            raise PropagationError(f"Sat {self.satid} {SGP4_ERRORS[status]}")
        rteme = np.array(position_eci)
        recef = np.empty(3, dtype=np.double)
        teme2ecef(mjd, rteme, recef)
        return recef


def unkozai(no_kozai, ecco, inclo, whichconst):
    """
    Undo Kozai transformation
    Ref: orbit_predictor/utils.py
    """
    _, _, _, xke, j2, _, _, _ = whichconst
    ak = pow(xke / no_kozai, 2.0 / 3.0)
    d1 = 0.75 * j2 * (3.0 * cos(inclo)**2 - 1.0) / (1.0 - ecco**2)**(3/2)
    del_ = d1 / ak ** 2
    adel = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
    return no_kozai / (1.0 + d1/adel**2)
