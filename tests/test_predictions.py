# Test using pytest

import pytest

from passpredict import (
    predict_single_satellite_overpasses,
)


# @pytest.mark.predict
# @pytest.mark.asyncio
# def test_predict_nocache():
#     """
#     Just confirm that the predict() function doesn't error
#     """
#     jd, location, sun_rECEF, satellite = init_find_overpasses
#     min_elevation = 10.0
#     location = Location(
#             name="",
#             latitude_deg=lat,
#             longitude_deg=lon,
#             elevation_m=h
#         )
#     # Get TLE data for satellite
#     tle = await tle_source.get_predictor(satid, today)

#     overpass_result = await run_in_threadpool(
#         predict_single_satellite_overpasses,
#         tle,
#         location,
#         date_start=today,
#         days=days,
#         min_elevation=10.0,
#     )
#     assert len(overpasses) > 0


# def test_predict_nocache_visible_only(init_predict):
#     """
#     Confirm that the predict() function only returns visible overpasses when requested
#     """
#     location, satellite, date_start, date_end = init_predict
#     overpasses = predictions.predict(location, satellite, date_start=date_start, date_end=date_end, dt_seconds=10, verbose=True, visible_only=True)
#     for overpass in overpasses:
#         assert overpass.type == PassType.visible


# Create a prediction test suite comparing results to
#    1. heavens-above.com
#    2. calsky.com
#    3. n2Y0.com
#    4. skyfield
#    5. pyephem
#
# Use multiple satellites and multiple locations at different latitudes
#   Satellites:
#     ISS (25544)
#     Hubble
#     Starlink-8
#     Lightsail 2
#     Envisat
#     Terra
#     a geosynchronus sat
#     a retrograde orbit sat
#     a molniya orbit sat
#     a sun-synchronus orbit sat
#     a geostationary orbit sat
#
#   Locations (lat, lon):
#     Longyearbyen, Norway (78.2)
#     Inuvik, Canada (68.3607 , -133.7230)
#     Helsinki, Finland (60.17)
#     Kiev, Ukraine (50.45)
#     New York, NY, USA (40.7128, -74.0060)
#     Austin, Texas, USA (30.2672, -97.7431)
#     Mexico City, Mexico (~20)
#     Bissau, Guinea Bissau (~10)
#     Quito, Ecuador (~0)
#     Johannesburg, South Africa (-26)
#     Sydney, Australia (-33)
#     Christchurch, New Zealand (-43)
#     Punta Arenas, Chile (-53)
#     McMurdo Station, Antarctica (-77)

# predict_results = {

# }