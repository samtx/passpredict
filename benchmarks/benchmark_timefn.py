# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.time import Time

from passpredict.predictions import compute_time_array
from passpredict.timefn import julian_date, jd2jc, jd2utc1, jday2datetime, invjday

class ComputeTimeArray:
    
    params = [1, 2, 5, 10, 14, 21]
    param_names = ['days']
    
    def setup(self, days):
        self.dt_start = datetime(2020, 7, 14, 11, 17, 00, tzinfo=timezone.utc)
        self.dt_end = self.dt_start + timedelta(days=days)
        

    def time_compute_time_array(self, *args):
        compute_time_array(self.dt_start, self.dt_end, 1.0)


class Jday:

    params = [2450383.09722222,  2453101.828154745, 2453101.8274067827]
    param_names = ['jd']

    def time_jday2datetime(self, jd):
        jday2datetime(jd)

    def time_invjday(self, jd):
        invjday(jd)