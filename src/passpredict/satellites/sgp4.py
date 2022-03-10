from __future__ import annotations
import typing
from functools import cached_property, lru_cache
from math import cos, radians, pi

import numpy as np
from sgp4.api import SGP4_ERRORS, Satrec
from sgp4.model import WGS84
from sgp4.earth_gravity import wgs84 as wgs84_params

from .base import SatellitePropagatorBase
from .._time import mjd2jdfr
from ..orbit import TLE, Orbit
from .._rotations import teme2ecef
from ..exceptions import PropagationError

if typing.TYPE_CHECKING:
    from ..sources import TLE


class SGP4Propagator(SatellitePropagatorBase):
    """
    SGP4 Propagator for satellite position and velocity.
    """
    def __init__(self, orbit: Orbit):
        """
        Params:
            satid: int = NORAD id for satellite
            source: PasspredictTLESource
        """
        self.orbit = orbit
        self.satid = orbit.satid
        self.name = orbit.name
        self.set_propagator(
            self.orbit.jdepoch,
            self.orbit.bstar,
            self.orbit.ndot,
            self.orbit.nddot,
            self.orbit.ecc,
            self.orbit.argp,
            self.orbit.inc,
            self.orbit.mo,
            self.orbit.no_kozai,
            self.orbit.raan
        )
        self.intrinsic_mag = 1.0   # ISS is -1.8



    @classmethod
    def from_tle(cls, tle: TLE) -> SGP4Propagator:
        orbit = Orbit.from_tle(tle)
        sat = cls(orbit)
        return sat

    def set_propagator(
        self,
        jdepoch: float,
        bstar: float,
        ndot: float,
        nddot: float,
        ecc: float,
        argp: float,
        inc: float,
        mo: float,
        no_kozai: float,
        raan: float,
    ):
        """
        sat.sgp4init(
            WGS72,           # gravity model
            'i',             # 'a' = old AFSPC mode, 'i' = improved mode
            5,               # satnum: Satellite number
            18441.785,       # epoch: days since 1949 December 31 00:00 UT
            2.8098e-05,      # bstar: drag coefficient (1/earth radii)
            6.969196665e-13, # ndot (NOT USED): ballistic coefficient (revs/day)
            0.0,             # nddot (NOT USED): mean motion 2nd derivative (revs/day^3)
            0.1859667,       # ecco: eccentricity
            5.7904160274885, # argpo: argument of perigee (radians)
            0.5980929187319, # inclo: inclination (radians)
            0.3373093125574, # mo: mean anomaly (radians)
            0.0472294454407, # no_kozai: mean motion (radians/minute)
            6.0863854713832, # nodeo: right ascension of ascending node (radians)
        )
        """
        xpdotp   =  1440.0 / (2.0 * pi);  #  229.1831180523293
        jd0 = 2433281.5   # julian date Dec 31 1949 00:00 UT
        satid = self.satid if isinstance(self.satid, int) else 0
        self._propagator = Satrec()
        self._propagator.sgp4init(
            WGS84,
            'i',
            satid,                  # satnum
            jdepoch - jd0,          # epoch
            bstar,                  # bstar
            ndot / (xpdotp * 1440.0),               # ndot
            nddot / (xpdotp * 1440.0 * 1440.0),     # nddot
            ecc,                    # ecco
            radians(argp),          # argpo
            radians(inc),           # inclo
            radians(mo),            # mo
            no_kozai / xpdotp,      # no_kozai
            radians(raan),          # nodeo in sgp4
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
