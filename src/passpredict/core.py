from __future__ import annotations
import datetime
import typing

from .observers import Observer
from .locations import Location

if typing.TYPE_CHECKING:
    from .satellites import SatellitePredictorBase


def predict_all_visible_satellite_overpasses(
    tles,
    location: Location,
    date_start: datetime.date,
    min_elevation: float,
):
    results = []
    for tle in tles:
        sat_results = predict_single_satellite_overpasses(
            tle, location,
            date_start=date_start,
            days=1,
            min_elevation=min_elevation
        )
        results.extend(sat_results)
    return results


def predict_single_satellite_overpasses(
    satellite: SatellitePredictorBase,
    location: Location,
    date_start: datetime.date,
    days: int,
    min_elevation: float
):
    date_end = date_start + datetime.timedelta(days=days)
    observer = Observer(location, satellite, aos_at_dg=min_elevation, tolerance_s=1.0)
    pass_iterator = observer.iter_passes(date_start, limit_date=date_end)
    passes = list(pass_iterator)
    return passes


def predict_next_overpass(
    satellite: SatellitePredictorBase,
    location: Location,
    date_start: datetime.date,
    min_elevation: float = 10.0
):
    observer = Observer(location, satellite, aos_at_dg=min_elevation, tolerance_s=1.0)
    pass_ = observer.get_next_pass(date_start)
    return pass_


def get_next_pass_detail(
    satellite: SatellitePredictorBase,
    location: Location,
    date_start: datetime.datetime,
    min_elevation: float = 10.0,
):
    observer = Observer(location, satellite, aos_at_dg=min_elevation, tolerance_s=1.0)
    pass_detail, llh = observer.get_next_pass_detail(date_start)
    return pass_detail, llh


def get_satellite_llh(
    satellite: SatellitePredictorBase,
    date_start: datetime.datetime,
    date_end: datetime.datetime,
    dt_seconds: float = 1,
):
    assert dt_seconds > 0, "dt_seconds must be greater than 0"
    time_step = datetime.timedelta(seconds=dt_seconds)
    n_steps = int((date_end - date_start).total_seconds() / dt_seconds) + 1
    llh = satellite.get_position_detail(date_start, n_steps, time_step)
    return llh