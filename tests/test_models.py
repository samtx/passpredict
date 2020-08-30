# test_models.py
import datetime

import numpy as np
from numpy.testing import assert_array_equal
import pytest
from astropy.time import Time

from passpredict import models
from passpredict.timefn import julian_date
from passpredict.propagate import compute_satellite_data
from passpredict.solar import compute_sun_data
from passpredict.schemas import Tle
from passpredict.tle import epoch_from_tle


# def test_find_overpasses_visible_only():
#     tle1 = "1 25544U 98067A   20154.57277630  .00016717  00000-0  10270-3 0  9118"
#     tle2 = "2 25544  51.6443  60.8122 0001995  12.6931 347.4269 15.49438452 29742"
#     tle = Tle(tle1=tle1, tle2=tle2, epoch=epoch_from_tle(tle1), satid=25544)
#     sat = compute_satellite_data(tle, )
#     store_sat_id = True if len(sats) > 1 else False
#     overpasses = []
#     for sat in sats:
#         rho = RhoVector(jd, sat, location, sun)
#         sat_overpasses = rho.find_overpasses(min_elevation, store_sat_id, visible_only=visible_only)
#         overpasses += sat_overpasses
#     return overpasses


# def test_Time_Model():
#     jd0 = julian_date(2020, 7, 1, 12, 34, 56)
#     jd = jd0 + np.linspace(0, 14, 100)
#     t = models.Time(jd=jd)
#     t2 = t[::2]
#     assert_array_equal(t.jd, jd)
#     assert_array_equal(t2.jd, jd[::2])


if __name__ == "__main__":
    pytest.main(['-v', __file__])