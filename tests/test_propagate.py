# import datetime
# import os

# import numpy as np
# import pytest

# from app import propagate
# from app.timefn import julian_date_array_from_datetime, jday2datetime_us_array
# from app.schemas import Satellite, Tle
# from app.utils import epoch_from_tle

# tz_utc = datetime.timezone.utc

# def compute_skyfield_ecef_position(tle1, tle2, jd):
#     """
#     Use the skyfield api to compute the position coordinates for satellite
#     """
#     from skyfield.api import Topos, load, EarthSatellite
#     from skyfield.sgp4lib import TEME_to_ITRF
#     ts = load.timescale()
#     sat = EarthSatellite(tle1, tle2, ts=ts)
#     datetimes = jday2datetime_us_array(jd)
#     t = ts.from_datetimes(datetimes)
#     rTEME, vTEME, _ = sat._position_and_velocity_TEME_km(t)
#     rECEF, _ = TEME_to_ITRF(jd, rTEME, vTEME)
#     return rECEF


# @pytest.mark.slow
# @pytest.mark.xfail
# def test_propagate_iss():
#     """
#     Compare results of propagate() with skyfield EarthSatellite.at()

#     datetime_start = 2020-06-01:00:00:00 UTC
#     datetime_end   = 2020-06-11:00:00:00 UTC
#     dt_seconds = 5

#     data/skyfield_iss_rECEF.npy
#     """
#     satellite = Satellite(id=25544, name="ISS")
#     tle1 = "1 25544U 98067A   20154.57277630  .00016717  00000-0  10270-3 0  9118"
#     tle2 = "2 25544  51.6443  60.8122 0001995  12.6931 347.4269 15.49438452 29742"
#     datetime_start = datetime.datetime(2020, 6, 1, 0, 0, 0, tzinfo=tz_utc)
#     datetime_end = datetime.datetime(2020, 6, 11, 0, 0, 0, tzinfo=tz_utc)
#     dt_sec = 30
#     tle = Tle.from_string(tle1, tle2)
#     jd = julian_date_array_from_datetime(datetime_start, datetime_end, dt_sec)
#     skyfield_rECEF = compute_skyfield_ecef_position(tle.tle1, tle.tle2, jd)
#     sat = propagate.compute_satellite_data(tle, jd)

#     # fname = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/skyfield_iss_rECEF.npy'))
#     # skyfield_rECEF = np.load(fname, allow_pickle=True)
#     diff = np.linalg.norm(sat.rECEF - skyfield_rECEF.T, axis=1)

#     # TODO: Make this more accurate in the future!!
#     np.testing.assert_array_less(diff, 40.0, verbose=True)  # difference less than 20 km


# if __name__ == "__main__":
#     import pytest
#     pytest.main(['-v', __file__])