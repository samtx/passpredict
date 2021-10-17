# # Test passpredict/schemas.py
# import datetime

# import pytest

# from app import schemas
# from app.utils import epoch_from_tle


# def test_Satellite():
#     satid = 25544
#     name = "International Space Station"
#     satellite = schemas.Satellite(id=satid, name=name)
#     assert satellite.id == satid
#     assert satellite.name == name


# # @pytest.mark.xfail(strict=True)
# def test_Satellite_fail():
#     """
#     This test shows that you must use keyword arguments for pydantic schema objects
#     """
#     satid = 25544
#     name = "International Space Station"
#     with pytest.raises(TypeError) as err:
#         satellite = schemas.Satellite(satid, name)
#     assert '__init__() takes' in str(err.value)


# def test_Point():
#     dt = datetime.datetime(2020, 6, 1, 12, 34, 56, tzinfo=datetime.timezone.utc)
#     timestamp = dt.timestamp()
#     azimuth = 270.5
#     elevation = 84.1
#     range_ = 3200.214
#     point = schemas.Point(
#         datetime=dt,
#         timestamp=timestamp,
#         azimuth=azimuth,
#         elevation=elevation,
#         range=range_
#     )
#     assert point.datetime == dt
#     assert point.timestamp == timestamp
#     assert point.azimuth == azimuth
#     assert point.elevation == elevation
#     assert point.range == range_


# def test_Overpass():
#     dt_start = datetime.datetime(2020, 6, 1, 12, 34, 56, tzinfo=datetime.timezone.utc)
#     start_pt = schemas.Point(
#         datetime = dt_start,
#         timestamp = dt_start.timestamp(),
#         azimuth = 270.5,
#         elevation = 12,
#         range = 3200.214
#     )
#     dt_max = datetime.datetime(2020, 6, 1, 12, 35, 59, tzinfo=datetime.timezone.utc)
#     max_pt = schemas.Point(
#         datetime = dt_max,
#         timestamp = dt_max.timestamp(),
#         azimuth = 358,
#         elevation = 84.1,
#         range = 2500
#     )
#     dt_end = datetime.datetime(2020, 6, 1, 12, 37, 12, tzinfo=datetime.timezone.utc)
#     end_pt = schemas.Point(
#         datetime = dt_end,
#         timestamp = dt_end.timestamp(),
#         azimuth = 20,
#         elevation = 10.0000001,
#         range = 3987.1483321548
#     )
#     overpass = schemas.Overpass(
#         start_pt=start_pt,
#         max_pt=max_pt,
#         end_pt=end_pt
#     )
#     assert overpass.start_pt == start_pt
#     assert overpass.start_pt.datetime == datetime.datetime(2020, 6, 1, 12, 34, 56, tzinfo=datetime.timezone.utc)
#     assert overpass.max_pt == max_pt


# def test_Location():
#     lat = 30.2672
#     lon = -97.7431
#     h = 88
#     name = 'Austin, Texas'
#     austin = schemas.Location(lat=lat, lon=lon, h=h, name=name)
#     assert austin.lat == lat
#     assert austin.lon == lon
#     assert austin.h == h
#     assert austin.name == name


# if __name__ == "__main__":
#     pytest.main(['-v', __file__])
