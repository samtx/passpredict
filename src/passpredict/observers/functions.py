from __future__ import annotations
import typing
import datetime

import numpy as np

from .core import PassType
from ..constants import R_EARTH
from ..time import julian_date

if typing.TYPE_CHECKING:
    from .observer import Observer


class VisualPoints(typing.NamedTuple):
    vis_begin_mjd: float = None
    vis_tca_mjd: float = None
    vis_end_mjd: float = None

if typing.TYPE_CHECKING:
    from .observer import Observer


class VisualPoints(typing.NamedTuple):
    vis_begin_mjd: float = None
    vis_tca_mjd: float = None
    vis_end_mjd: float = None


def make_utc(d: datetime.datetime) -> datetime.datetime:
    """ Make a datetime a UTC timezone aware datetime """
    if not d:
        return d
    if d.tzinfo:
        res = d.astimezone(datetime.timezone.utc)
    else:
        res = d.replace(tzinfo=datetime.timezone.utc)
    return res


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


def visual_pass_details(
    observer: Observer,
    aos_mjd: float,
    tca_mjd: float,
    los_mjd: float,
    *,
    tol: float = 1/86400.0,
    sunrise_dg: float = -6,
    n: int = 5,  # number of steps
) -> typing.Tuple([PassType, VisualPoints]):
    t0 = aos_mjd
    tf = los_mjd
    t = np.linspace(t0, tf, n)

    def sun_el_fn(t):
        return observer.location._sun_elevation_mjd(t) - sunrise_dg

    el = np.array([sun_el_fn(t_) for t_ in t])

    if np.min(el) > 0:
        # entire pass in sunlit
        return (PassType.daylight, VisualPoints(None, None, None))

    if el[0]*el[-1] < 0:
        # part of the pass is in darkness.
        # only part of the pass is sunlit. Find new jd0, jdf
        x = find_root(sun_el_fn, t0, tf, tol)
        tmp1 = sun_el_fn(x - tol)
        tmp2 = sun_el_fn(x + tol)
        if tmp1 < tmp2:
            # sun elevation is decreasing
            t0 = x
        else:
            tf = x
    # Now use t0 and tf to find when satellite is
    # illuminated by sun
    t = np.linspace(t0, tf, n)

    def illum_fn(t):
        return observer.satellite._illumination_distance_mjd(t) - R_EARTH  # noqa

    illum_pts = np.array([illum_fn(t_) for t_ in t])

    if np.max(illum_pts) < 0:
        # entire pass is in shadow
        return (PassType.unlit, VisualPoints(None, None, None))

    if illum_pts[0]*illum_pts[-1] < 0:
        # the satellite is visible for only part of the pass.
        # Find new t0, tf
        x = find_root(illum_fn, t0, tf, tol)
        tmp1 = illum_fn(x - tol)
        tmp2 = illum_fn(x + tol)
        if tmp1 < tmp2:
            # satellite is coming out of shadow
            t0 = x
        else:
            # satellite is going into shadow
            tf = x

    # Set visible start and end points for Pass
    vis_begin_mjd = t0
    vis_end_mjd = tf
    # Find maximum elevation during visible period
    if vis_begin_mjd <= tca_mjd <= vis_end_mjd:
        vis_tca_mjd = tca_mjd
    elif observer._elevation_mjd(vis_begin_mjd) > observer._elevation_mjd(vis_end_mjd):
        vis_tca_mjd = vis_begin_mjd
    else:
        vis_tca_mjd = vis_end_mjd
    visual_points = VisualPoints(vis_begin_mjd, vis_tca_mjd, vis_end_mjd)
    return (PassType.visible, visual_points)
