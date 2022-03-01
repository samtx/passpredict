from __future__ import annotations
import datetime
import typing
from math import radians, pi

import numpy as np
from scipy.optimize import root_scalar

from .core import BasicPassInfo, PassType
from .functions import make_utc
from .._time import datetime2mjd
from ..constants import R_EARTH
from ..exceptions import PropagationError

if typing.TYPE_CHECKING:
    from .observer import Observer


def orbit_predictor_iterator(
    observer: Observer,
    start_date: datetime.datetime,
    limit_date: datetime.datetime = None,
    *,
    aos_at_dg: float = 0,
    max_elevation_gt: float = 0,
    tol: float = 1,  # tolerance in seconds
    sunrise_dg: float = -6,
    **extra_options,
):
    """Returns one pass each time"""
    start_date = make_utc(start_date)
    limit_date = make_utc(limit_date)

    # Get extra options
    aos_at = radians(aos_at_dg)
    max_elevation_gt = radians(max([max_elevation_gt, aos_at_dg]))
    if tol <= 0:
        raise Exception("Tolerance must be > 0")
    tol = tol/86400.0
    mjd = datetime2mjd(start_date)
    limit_mjd = datetime2mjd(limit_date) if limit_date else None

    while True:
        if limit_mjd is not None and mjd > limit_mjd:
            break

        if _is_ascending(observer, mjd, tol):
            # we need a descending point
            ascending_mjd = mjd

            # descending_date = _find_nearest_descending(observer, ascending_date)
            # def _find_nearest_descending(observer, ascending_date, tolerance):
            candidate_found = False
            for candidate in _sample_points(observer, ascending_mjd):
                if not _is_ascending(observer, candidate, tol):
                    descending_mjd = candidate
                    candidate_found = True
                    break
            if not candidate_found:
                # logger.error('Could not find a descending pass over %s start date: %s - TLE: %s',
                #              self.location, ascending_date, self.satellite.tle)
                raise Exception("Can not find an descending phase")

            # Find TCA tca_dt = _find_tca(observer, ascending_date, descending_date)
            while not (descending_mjd - ascending_mjd <= tol):  # precision reached
                midpoint = _midpoint(ascending_mjd, descending_mjd)
                if _is_ascending(observer, midpoint, tol):
                    ascending_mjd = midpoint
                else:
                    descending_mjd = midpoint
            tca_mjd = ascending_mjd

            tca_elevation = observer._elevation_at_mjd(tca_mjd)
            if tca_elevation > max_elevation_gt:

                # Find AOS
                end = tca_mjd
                start = tca_mjd - _orbit_step(observer, 0.34)  # On third of the orbit
                elevation = observer._elevation_at_mjd(start)
                while not (end - start <= tol):  # precision reached
                    midpoint = _midpoint(start, end)
                    elevation = observer._elevation_at_mjd(midpoint)
                    if elevation < aos_at:
                        start = midpoint
                    else:
                        end = midpoint
                aos_mjd = end

                # Find LOS  los_dt = self._find_los(tca_dt)
                start = tca_mjd
                end = tca_mjd + _orbit_step(observer, 0.34)
                while not (end - start <= tol):  # precision reached
                    midpoint = _midpoint(start, end)
                    elevation = observer._elevation_at_mjd(midpoint)
                    if elevation < aos_at:
                        end = midpoint
                    else:
                        start = midpoint
                los_mjd = start

            else:
                mjd = tca_mjd + _orbit_step(observer, 0.6)
                continue

            # Find visual pass details
            # First, get endpoints of when location is not sunlit
            t0 = aos_mjd
            tf = los_mjd
            t = np.linspace(t0, tf, 5)  # use 5 points for spline

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
                    result = root_scalar(
                        sun_el_fn, method='bisect', bracket=(t0, tf),
                        x0=t0, xtol=tol
                    )
                    x = result.root
                    tmp1 = sun_el_fn(x - tol)
                    tmp2 = sun_el_fn(x + tol)
                    if tmp1 < tmp2:
                        # sun elevation is decreasing
                        t0 = x
                    else:
                        tf = x

                # Now use jd0 and jdf to find when satellite
                # is illuminated by sun
                t = np.linspace(t0, tf, 5)  # use 5 points for spline

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
                        result = root_scalar(
                            illum_fn, method='bisect', bracket=(t0, tf),
                            x0=t0, xtol=tol
                        )
                        x = result.root
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
            mjd = pass_.tca + _orbit_step(observer, 0.6)

            if limit_mjd is not None and pass_.aos > limit_mjd:
                break
        else:
            candidate_found = False
            for candidate in _sample_points(observer, mjd):
                if _is_ascending(observer, candidate, tol):
                    mjd = candidate
                    candidate_found = True
                    break
            if not candidate_found:
                msg = (
                    f'Sat {observer.satellite.satid},'
                    f'date: {mjd}, could not find an ascending phase'
                )
                raise Exception(msg)


def _sample_points(observer: Observer, mjd: float):
    """Helper method to found ascending or descending phases of elevation"""
    start = mjd
    end = mjd + _orbit_step(observer, 0.99)
    mid = _midpoint(start, end)
    mid_right = _midpoint(mid, end)
    mid_left = _midpoint(start, mid)
    mid_right_2 = _midpoint(mid, mid_right)
    mid_left_2 = _midpoint(start, mid_left)
    mid_right_3 = _midpoint(mid_right, end)
    mid_left_3 = _midpoint(mid_left, mid)
    pts = [
        end, mid, mid_right, mid_left, mid_right_2,
        mid_left_2, mid_right_3, mid_left_3
    ]
    return pts


def _midpoint(start, end):
    """Returns the midpoint between two dates"""
    return start + (end - start) / 2


def _is_ascending(observer: Observer, mjd: float, tol: float):
    """Check is elevation is ascending or descending on a given point"""
    elevation = observer._elevation_at_mjd(mjd)
    next_elevation = observer._elevation_at_mjd(mjd + tol)
    return elevation <= next_elevation


def _orbit_step(observer: Observer, size: float) -> float:
    """
    Returns a mjd time step, that will make the satellite
    advance a given number of orbits
    """
    step_in_radians = size * 2 * pi
    seconds = (step_in_radians / observer.satellite.mean_motion) * 60
    return seconds / 86400.0
