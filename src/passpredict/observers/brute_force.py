from __future__ import annotations
import datetime
from typing import TYPE_CHECKING
from math import radians, degrees

from .core import BasicPassInfo
from .functions import make_utc, find_root, find_min, visual_pass_details
from .._time import datetime2mjd

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
                return observer._elevation_mjd(t) - aos_at

            aos_mjd = find_root(el_fn, prev_mjd, mjd, tol)
            mjd = aos_mjd + step

            # find los
            while observer._elevation_mjd(mjd) > aos_at:
                prev_mjd = mjd
                mjd += step

            los_mjd = find_root(el_fn, prev_mjd, mjd, tol)

            # find tca
            def el_fn(t):
                return -observer._elevation_mjd(t)

            tca_mjd, max_el_rad = find_min(el_fn, aos_mjd, los_mjd, tol)

            tca_elevation = degrees(-max_el_rad)

            # Find visual pass details
            type_, visual_points = visual_pass_details(
                observer,
                aos_mjd,
                tca_mjd,
                los_mjd,
                tol=tol,
                sunrise_dg=sunrise_dg,
                n=25,
            )
            pass_ = BasicPassInfo(
                aos_mjd,
                tca_mjd,
                los_mjd,
                tca_elevation,
                type_=type_,
                vis_begin_mjd=visual_points.vis_begin_mjd,
                vis_end_mjd=visual_points.vis_end_mjd,
                vis_tca_mjd=visual_points.vis_tca_mjd,
            )
            yield pass_
            mjd = pass_.los_mjd + step*5
        if limit_mjd and mjd > limit_mjd:
            break
        prev_mjd = mjd
        mjd += step


def _above_horizon(observer, t, aos_at):
    el = observer._elevation_mjd(t)
    return el >= aos_at


def _crosses_horizon(observer, t1, t2, aos_at):
    if not _above_horizon(observer, t1, aos_at) and _above_horizon(observer, t2, aos_at):
        return True
    else:
        return False
