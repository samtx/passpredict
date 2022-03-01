from __future__ import annotations
import datetime
from telnetlib import EC
from typing import TYPE_CHECKING
from math import radians, degrees

import numpy as np

from .core import BasicPassInfo, PassType
from .functions import make_utc, find_root, find_min
from .._time import datetime2mjd
from ..constants import R_EARTH


if TYPE_CHECKING:
    from .observer import Observer


def brute_force_iterator(
    observer: Observer,
    start_date: datetime.datetime,
    limit_date: datetime.datetime = None,
    *,
    aos_at_dg: float = 0,
    max_elevation_gt: float = 0,
    tol: float = 1,  # tolerance in seconds
    sunrise_dg: float = -6,
    time_step: float = 10,  # time step in seconds
    **extra_options,
) -> BasicPassInfo:
    """
    Predicts passes of a satellite over a given location.
    Exposes an iterable interface.
    This is a brute force observation algorithm, useful for validation

    Returns one pass each time
    """
    start_date = make_utc(start_date)
    limit_date = make_utc(limit_date)
    mjd = datetime2mjd(start_date)
    limit_mjd = datetime2mjd(limit_date) if limit_date else None

    # Get extra options
    aos_at = radians(aos_at_dg)
    max_elevation_gt = radians(max([max_elevation_gt, aos_at_dg]))
    if tol <= 0:
        raise Exception("Tolerance must be > 0")
    if tol >= time_step:
        raise Exception("Tolerance must be greater than time_step")
    tol = tol / 86400.0
    step = time_step / 86400.0

    aos_at = radians(aos_at_dg)

    prev_mjd = mjd - step
    while True:
        if _crosses_horizon(observer, prev_mjd, mjd, aos_at):
            # satellite has just come above the horizon, find aos, tca, and los

            def el_fn(t):
                return observer._elevation_at_mjd(t) - aos_at

            aos_mjd = find_root(el_fn, prev_mjd, mjd, tol)
            mjd = aos_mjd + step

            # find los
            while observer._elevation_at_mjd(mjd) > aos_at:
                prev_mjd = mjd
                mjd += step

            los_mjd = find_root(el_fn, prev_mjd, mjd, tol)

            # find tca
            def el_fn(t):
                return -observer._elevation_at_mjd(t)

            tca_mjd, max_el_rad = find_min(el_fn, aos_mjd, los_mjd, tol)

            tca_elevation = degrees(-max_el_rad)

            # Find visual pass details
            # First, get endpoints of when location is not sunlit
            t0 = aos_mjd
            tf = los_mjd
            t = np.linspace(t0, tf, 25)

            def sun_el_fn(t):
                return observer.location.sun_elevation_mjd(t) - sunrise_dg

            el = np.array([sun_el_fn(t_) for t_ in t])

            pass_ = None
            if np.min(el) > 0:
                # entire pass in sunlit
                pass_ = BasicPassInfo(
                    aos_mjd, tca_mjd, los_mjd, tca_elevation,
                    type_=PassType.daylight
                )

            else:
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
                t = np.linspace(t0, tf, 25)

                def illum_fn(t):
                    return observer.satellite.illumination_distance_mjd(t) - R_EARTH  # noqa

                illum_pts = np.array([illum_fn(t_) for t_ in t])

                if np.max(illum_pts) < 0:
                    # entire pass is in shadow
                    pass_ = BasicPassInfo(
                        aos_mjd, tca_mjd, los_mjd, tca_elevation,
                        type_=PassType.unlit
                    )

                else:
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
                    elif observer._elevation_at_mjd(vis_begin_mjd) > observer._elevation_at_mjd(vis_end_mjd):
                        vis_tca_mjd = vis_begin_mjd
                    else:
                        vis_tca_mjd = vis_end_mjd
                    pass_ = BasicPassInfo(
                        aos_mjd, tca_mjd, los_mjd, tca_elevation,
                        type_=PassType.visible, vis_begin_mjd=vis_begin_mjd,
                        vis_end_mjd=vis_end_mjd, vis_tca_mjd=vis_tca_mjd,
                    )
            yield pass_
            mjd = pass_.los_mjd + step*5
        if limit_mjd and mjd > limit_mjd:
            break
        prev_mjd = mjd
        mjd += step


def _above_horizon(observer, t, aos_at):
    el = observer._elevation_at_mjd(t)
    return el >= aos_at


def _crosses_horizon(observer, t1, t2, aos_at):
    if not _above_horizon(observer, t1, aos_at) and _above_horizon(observer, t2, aos_at):
        return True
    else:
        return False
