from __future__ import annotations
import datetime
import time
import typing
from math import radians, degrees, pi, log10, sin, cos, acos, floor
from enum import Enum

import numpy as np
# from scipy.optimize import minimize_scalar

from .base import ObserverBase, BasicPassInfo
from .functions import _make_utc, julian_date_sum
from ..time import julian_date, julian_date_from_datetime
from .._time import jday2datetime_us
from ..exceptions import NotReachable, PropagationError
from ..locations import Location

if typing.TYPE_CHECKING:
    from ..satellites import SatellitePredictor, LLH


class BruteForceObserver(ObserverBase):
    """
    Predicts passes of a satellite over a given location.
    Exposes an iterable interface.
    This is a brute force observation algorithm, useful for validation
    """

    def __init__(
        self,
        location: Location,
        satellite: SatellitePredictor,
        time_step=10,
        **kwargs
    ):
        """
        Initialize Observer but also compute radians for geodetic coordinates
        """
        super().__init__(location, satellite, **kwargs)
        if time_step <= 0:
            raise Exception("Time step must be > 0")
        self.jd_step = time_step / 86400

    def iter_passes(self, start_date, limit_date=None, visible_only=False):
        """Returns one pass each time"""
        start_date = _make_utc(start_date)
        limit_date = _make_utc(limit_date)
        jd = julian_date_sum(start_date)
        if not limit_date:
            limit_jd = None
        else:
            limit_jd = julian_date_sum(limit_date)
        # First check if satellite is currently above horizon
        # If so, rewind start time until below horizon
        # stop rewinding after 1 day
        prev_jd_limit = jd - 1
        while self._above_horizon(jd) and jd > prev_jd_limit:
            jd = jd - self.jd_step

        prev_jd = jd - self.jd_step
        while True:
            if self._crosses_horizon(prev_jd, jd):
                # satellite has just come above the horizon, find aos, tca, and los
                pass_, los_jd = self._refine_pass(prev_jd, jd)
                if self._is_pass_valid(pass_, visible_only=visible_only):
                    predicted_pass = self._build_predicted_pass(pass_)
                    yield predicted_pass
                    jd = los_jd + self.jd_step * 5
            if limit_jd and jd > limit_jd:
                break
            prev_jd = jd
            jd += self.jd_step

    def set_minimum_elevation(self, elevation: float):
        """  Set minimum elevation for an overpass  """
        self.aos_at = radians(elevation)
        self.aos_at_deg = elevation

    def _above_horizon(self, jd):
        el = self._elevation_at_jd(jd)
        return el >= self.aos_at

    def _crosses_horizon(self, jd1, jd2):
        if not self._above_horizon(jd1) and self._above_horizon(jd2):
            return True
        else:
            return False

    def _refine_pass(self, jd1, jd2) -> BasicPassInfo:
        el_fn = lambda x: self._elevation_at_jd(x) - self.aos_at
        aos_jd = find_root(el_fn, jd1, jd2, self.jd_tol)
        jd = aos_jd + self.jd_step
        # find los
        while self._elevation_at_jd(jd) > self.aos_at:
            prev_jd = jd
            jd += self.jd_step
        los_jd = find_root(el_fn, prev_jd, jd, self.jd_tol)
        # find tca
        el_fn = lambda x: -self._elevation_at_jd(x)
        tca_jd, max_el_rad = find_min(el_fn, aos_jd, los_jd, self.jd_tol)
        max_el = degrees(-max_el_rad)
        aos_dt = jday2datetime_us(aos_jd)
        tca_dt = jday2datetime_us(tca_jd)
        los_dt = jday2datetime_us(los_jd)
        return BasicPassInfo(aos_dt, tca_dt, los_dt, max_el), los_jd


def find_root(f: typing.Callable, a: float, b: float, tol: float) -> float:
    """
    Find root of function f() with bisection method within tolerance tol.
    Return root.
    """
    assert a < b
    fa, fb = f(a), f(b)
    if fa*fb >= 0:
        return None
        # # return value closest to zero
        # if abs(fa) < abs(fb):
        #     return a
        # else:
        #     return b

    diff = tol + 1
    while diff > tol:
        mid = (a + b) / 2
        fmid = f(mid)
        if fa*fmid < 0:
            b = mid
            fb = f(b)
        elif fb*fmid < 0:
            a = mid
            fa = f(a)
        diff = abs(b - a)
    return mid


def find_min(f: typing.Callable, a: float, b: float, tol: float) -> float:
    """
    Find minimum of bounded univariate scalar function using scipy.optimize.minimize_scalar
    """
    assert a < b
    # res = minimize_scalar(f, bounds=(a, b), method='golden', tol=tol, options={'xatol': tol})
    diff = b - a
    fvec = np.vectorize(f)
    N = 5
    while diff > tol:
        x = np.linspace(a, b, N)
        y = fvec(x)
        i = y.argmin()
        if i == N:
            a = x[N - 1]
            b = x[N]
        elif i == 0:
            a = x[0]
            b = x[1]
        else:
            a = x[i - 1]
            b = x[i + 1]
        diff = abs(b - a)
    xsol = a + diff / 2
    fsol = f(xsol)
    return xsol, fsol



def julian_date_sum(d: datetime.datetime) -> float:
    jd, jdfr = julian_date(d.year, d.month, d.day, d.hour, d.minute, d.second + d.microsecond/1e6)
    return jd + jdfr
