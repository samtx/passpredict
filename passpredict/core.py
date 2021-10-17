import datetime

from .predictors import SatellitePredictor
from .locations import Location


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
    predictor: SatellitePredictor,
    location: Location,
    date_start: datetime.date,
    days: int,
    min_elevation: float
):
    date_end = date_start + datetime.timedelta(days=days)
    pass_iterator = predictor.pass_iterator(
        location, when_utc=date_start, limit_date=date_end,
        aos_at_dg=min_elevation, tolerance_s=1.0
    )
    passes = list(pass_iterator)
    return passes


def predict_next_overpass(
    predictor: SatellitePredictor,
    location: Location,
    date_start: datetime.date,
    min_elevation: float = 10.0
):
    pass_ = predictor.get_next_pass(
        location, aos_dt=date_start,
        aos_at_dg=min_elevation
    )
    return pass_