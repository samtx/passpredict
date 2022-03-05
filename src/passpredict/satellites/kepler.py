from __future__ import annotations
import datetime
from math import radians, sqrt, pi, sin, cos, asin, acos
from re import U
import typing
from functools import cached_property

import numpy as np
from orbit_predictor.keplerian import coe2rv as _coe2rv
from orbit_predictor.predictors.numerical import pkepler as _pkepler
from orbit_predictor.utils import mean_motion

from .base import SatellitePropagatorBase
from .._rotations import teme2ecef
from ..orbit import Orbit, TLE, kepler_eqn_E, nu_to_anomaly, anomaly_to_nu
from ..constants import MU_E, R_EARTH, J2_E, MJD0


def coe2rv(mu, p, ecc, inc, raan, argp, nu):
    return _coe2rv(mu, p, ecc, inc, raan, argp, nu)


def pkepler(dt, argp0, e0, inc, raan0, a0, nu0, ndot=0, nddot=0):
    p0 = a0 * (1-e0*e0)
    r, v = _pkepler(argp0, dt, e0, inc, p0, raan0, a0, nu0)
    return r, v
    # # Initial mean anomaly
    # if e0 != 0:
    #     E0 = nu_to_anomaly(nu0, e0)
    # else:
    #     # E0 = u
    #     raise Exception('eccentricity == 0')
    # M0 = E0 - e0 * sin(E0)
    # p0 = a0 * (1-e0*e0)
    # n0 = sqrt(MU_E / (a0*a0*a0))
    # # k = 3 * n0 * R_EARTH * R_EARTH * J2_E / (2 * p0 * p0)
    # # update for perturbation
    # a = a0 - (2*a0)/(3*n0) * ndot * dt
    # e = e0 - 2*(1 - e0)/(3*n0) * ndot * dt
    # raan = raan0 - 3*n0*R_EARTH*R_EARTH*J2_E/(2*p0*p0) * cos(inc) * dt
    # argp = argp0 + n0*R_EARTH*R_EARTH*J2_E/(4*p0*p0) * (4 - 5*sin(inc)*sin(inc)) * dt
    # M = M0 + n0*dt + ndot/2*dt*dt + nddot/6*dt*dt*dt
    # p = a * (1 - e*e)
    # E = kepler_eqn_E(M, e)
    # if e != 0:
    #     nu = anomaly_to_nu(E, e)
    # else:
    #     raise Exception('eccentricity == 0')
    # r, v = coe2rv(MU_E, p, e, inc, raan, argp, nu)
#     return r, v


def pkepler2(dt, argp0, e0, inc, raan0, a0, nu0, ndot=0, nddot=0):
    # Initial mean anomaly
    if e0 != 0:
        E0 = nu_to_anomaly(nu0, e0)
    else:
        # E0 = u
        raise Exception('eccentricity == 0')
    M0 = E0 - e0 * sin(E0)
    p0 = a0 * (1-e0*e0)
    n0 = sqrt(MU_E / (a0*a0*a0))
    # k = 3 * n0 * R_EARTH * R_EARTH * J2_E / (2 * p0 * p0)
    # update for perturbation
    a = a0 - (2*a0)/(3*n0) * ndot * dt
    e = e0 - 2*(1 - e0)/(3*n0) * ndot * dt
    raan = raan0 - 3*n0*R_EARTH*R_EARTH*J2_E/(2*p0*p0) * cos(inc) * dt
    argp = argp0 + n0*R_EARTH*R_EARTH*J2_E/(4*p0*p0) * (4 - 5*sin(inc)*sin(inc)) * dt
    M = M0 + n0*dt + ndot/2*dt*dt + nddot/6*dt*dt*dt
    p = a * (1 - e*e)
    E = kepler_eqn_E(M, e)
    if e != 0:
        nu = anomaly_to_nu(E, e)
    else:
        raise Exception('eccentricity == 0')
    r, v = coe2rv(MU_E, p, e, inc, raan, argp, nu)
    return r, v


class KeplerPropagator(SatellitePropagatorBase):
    """
    Propagator that uses the Keplerian osculating orbital elements.
    We use a naïve propagation algorithm that advances the anomaly the
    corresponding amount depending on the time difference and keeps all
    the rest of the osculating elements. It's robust against singularities
    as long as the starting elements are well specified but only works for
    elliptical orbits (ecc < 1). This limitation is not a problem since the object
    of study are artificial satellites orbiting the Earth.
    """
    def __init__(self, orbit: Orbit):
        """Initializes predictor.
        :param sma: Semimajor axis, km
        :param ecc: Eccentricity
        :param inc: Inclination, deg
        :param raan: Right ascension of the ascending node, deg
        :param argp: Argument of perigee, deg
        :param ta: True anomaly, deg
        :param epoch: Epoch, datetime
        """
        self.orbit = orbit
        self.satid = orbit.satid
        self.name = orbit.name
        epoch_jd = orbit.jdepoch + orbit.jdepochF
        self.epoch_mjd = epoch_jd - MJD0

    @property
    def sate_id(self):
        # Keplerian predictors are not made of actual observations
        return "<custom>"

    @property
    def mean_motion(self):
        """Mean motion, in radians per minute"""
        return mean_motion(self._sma) * 60

    @classmethod
    def from_tle(cls, tle: TLE):
        """Returns approximate keplerian elements from TLE.
        The conversion between mean elements in the TEME reference
        frame to osculating elements in any standard reference frame
        is not well defined in literature (see Vallado 3rd edition, pp 236 to 240)
        """
        orbit = Orbit.from_tle(tle)
        sat = cls(orbit)
        return sat

    @classmethod
    def from_coe(
        cls,
        n,
        ecc,
        inc,
        argp,
        raan,
        M,
        jdepoch,
        ndot=0,
        nddot=0
    ):
        orbit = Orbit.from_coe(n, ecc, inc, argp, raan, M, jdepoch, ndot, nddot)
        sat = cls(orbit)
        return sat

    @classmethod
    def from_rv(
        cls,
        r,
        v,
        jdepoch,
        ndot=0,
        nddot=0
    ):
        orbit = Orbit.from_rv(r, v, jdepoch, ndot, nddot)
        sat = cls(orbit)
        return sat

    @cached_property
    def coe_init(self):
        return (
            radians(self.orbit.argp),
            self.orbit.ecc,
            radians(self.orbit.inc),
            radians(self.orbit.raan),
            self.orbit.a,
            radians(self.orbit.nu),
            self.orbit.ndot,
            self.orbit.nddot,
        )

    def _position_ecef_mjd(self, mjd: float) -> np.ndarray:
        """ Use modified Julian date """

        dt = (mjd - self.epoch_mjd) * 86400.0  # dt in seconds
        argp0, e0, inc, raan0, a0, nu0, ndot, nddot = self.coe_init
        reci, _ = pkepler(dt, argp0, e0, inc, raan0, a0, nu0, ndot=ndot, nddot=nddot)
        recef = np.empty(3, dtype=np.double)
        teme2ecef(mjd, reci, recef)
        return recef


    def propagate_eci(self, when_utc=None):
        """Return position and velocity in the given date using ECI coordinate system.
        """
        if when_utc is None:
            when_utc = datetime.datetime.now(tz=datetime.timezone.utc)

        # Orbit parameters
        sma = self._sma
        ecc = self._ecc
        p = sma * (1 - ecc ** 2)
        inc = radians(self._inc)
        raan = radians(self._raan)
        argp = radians(self._argp)
        ta = radians(self._ta)

        delta_t_sec = (when_utc - self._epoch).total_seconds()

        # Propagate
        position_eci, velocity_eci = kepler(dt, argp, ecc, inc, p, raan, sma, ta)

        return tuple(position_eci), tuple(velocity_eci)