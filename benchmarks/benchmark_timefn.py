# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from datetime import datetime, timedelta, timezone

import numpy as np

from passpredict import time as pptime
from passpredict import _time

class Jday2Datetime:
    params = [2450383.09722222,  2453101.828154745, 2453101.8274067827]
    param_names = ['jd']

    def setup(self, jd):
        self.jd = 2453101.828154745

    def time_jday2datetime(self, jd):
        _time.jday2datetime(self.jd)
