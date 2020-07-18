# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from datetime import datetime, timezone, timedelta

import numpy as np
from astropy.time import Time

from passpredict.predictions import compute_sun_data, compute_time_array
from passpredict.models import SpaceObject, RhoVector, Sun
from passpredict.timefn import julian_date, jd2jc, jd2utc1, jday2datetime
from passpredict.solar import sun_pos, is_sat_illuminated


class SunPosition:
    
    params = [1, 2, 5, 10, 14, 21]
    param_names = ['days']
    
    def setup(self, days):
        dt_start = datetime(2020, 7, 14, 11, 17, 00, tzinfo=timezone.utc)
        dt_end = dt_start + timedelta(days=days)
        self.t = compute_time_array(dt_start, dt_end, 1.0)

    def time_compute_sun_data(self, *args):
        sun = compute_sun_data(self.t)


    def peakmem_compute_sun_data(self, *args):
        sun = compute_sun_data(self.t)
