# from datetime import datetime, timedelta, timezone
# import json
# import pathlib

# import pytest
# from pytest import approx

# from passpredict import Observer, Location, SGP4Predictor
# from passpredict import *
# from passpredict.observers import PassPoint, Visibility

# try:
#     from zoneinfo import ZoneInfo
# except ImportError:
#     from backports.zoneinfo import ZoneInfo




# def assert_datetime_approx(dt1, dt2, delta_seconds):
#     """
#     Compare two python datetimes. Assert difference is <= delta_seconds
#     """
#     diff = (dt1 - dt2).total_seconds()
#     assert diff == approx(0.0, abs=delta_seconds)


# def assert_passpoint_approx(pt1, pt2, *, dt_tol=1, el_tol=1, az_tol=10, range_tol=10, bright_tol=3):
#     """
#     Compare two PassPoint objects
#     """
#     assert_datetime_approx(pt1.dt, pt2.dt, dt_tol)
#     assert pt1.elevation == approx(pt2.elevation, abs=el_tol)
#     assert pt1.azimuth == approx(pt2.azimuth, abs=az_tol)
#     assert pt1.direction == pt2.direction
#     # if (pt1.brightness is not None) or (pt2.brightness is not None):
#     #     assert pt1.brightness == approx(pt2.brightness, abs=bright_tol)


# def test_n2yo_zurich_iss_visibility_predictions():
#     """
#     From n2yo.com. Queried 1/28/2022
#     """
#     data_path = pathlib.Path(__file__).parent / 'data' / 'n2yo_zurich_iss.json'
#     with open(data_path, encoding="utf-8") as f:
#         data = json.load(f)

#     tle_lines = tuple(data['tle'].splitlines())
#     tle = TLE(25544, tle_lines, name='ISS')
#     satellite = SGP4Predictor.from_tle(tle)
#     satellite.intrinsic_mag = -2.5  # From McCants qs.mag
#     location = Location("Zurich", 47.3744, 8.5410, 0)
#     observer = Observer(location, satellite, aos_at_dg=0, sunrise_dg=0, tolerance_s=0.25)


#     TZ = ZoneInfo('Europe/Zurich')
#     UTC = timezone.utc
#     start = datetime(2022, 1, 29, 0, 0, 0, tzinfo=UTC)
#     end = datetime(2022, 2, 6, tzinfo=UTC)
#     date = start
#     for n2yo_pass in data['passes']:
#         aos = PassPoint(
#             dt=datetime.fromtimestamp(n2yo_pass["startUTC"], tz=UTC),
#             elevation=n2yo_pass['startEl'],
#             azimuth=n2yo_pass['startAz'],
#             range=None,  # range is not provided as part of N2YO api
#         )
#         tca = PassPoint(
#             dt=datetime.fromtimestamp(n2yo_pass["maxUTC"], tz=UTC),
#             elevation=n2yo_pass['maxEl'],
#             azimuth=n2yo_pass['maxAz'],
#             range=None,  # range is not provided as part of N2YO api
#         )
#         los = PassPoint(
#             dt=datetime.fromtimestamp(n2yo_pass["endUTC"], tz=UTC),
#             elevation=n2yo_pass['endEl'],
#             azimuth=n2yo_pass['endAz'],
#             range=None,  # range is not provided as part of N2YO api
#         )

#         # Find next pass
#         tols = {
#             'dt_tol': 60,
#             'el_tol': 5,
#             'az_tol': 10,
#             'range_tol': 10,
#             'bright_tol': 3,
#         }
#         pass_ = observer.get_next_pass(date, limit_date=end, visible_only=True)
#         assert_passpoint_approx(pass_.vis_begin, aos, **tols)
#         assert_passpoint_approx(pass_.tca, tca, **tols)
#         assert_passpoint_approx(pass_.vis_end, los, **tols)
#         assert isinstance(pass_.brightness, float)
#         assert pass_.type == Visibility.visible
#         date = pass_.los.dt + timedelta(minutes=10)


# if __name__ == "__main__":
#     pytest.main([__file__])