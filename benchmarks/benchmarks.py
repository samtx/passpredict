# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import datetime

import numpy as np

import passpredict
from passpredict.schemas import Location, Satellite
from passpredict.predictions import predict_passes
from passpredict.propagate import propagate
from passpredict.timefn import truncate_datetime
from passpredict.utils import get_TLE

class Suite:
    """
    An example benchmark that times the performance of various kinds
    of iterating over dictionaries in Python.
    """
    def setup(self):
        # Set up satellite position
        dt_seconds = 1
        num_days = 14
        min_elevation = 10.0
        austin = Location(lat=30.2672, lon=-97.7431, h=0, name='Austin')
        iss = Satellite(id=25544, name='ISS')
        iss_tle = get_TLE(iss)
        datetime_start = truncate_datetime(
            datetime.datetime.now(tz=datetime.timezone.utc)
        )
        datetime_end = datetime_start + datetime.timedelta(days=num_days)
        iss_rv = propagate.__wrapped__(
            iss_tle.tle1, iss_tle.tle2, datetime_start, datetime_end, dt_seconds
        )
        self.austin = austin
        self.iss_rv = iss_rv
        self.min_elevation = min_elevation
        
    def time_predict_passes(self):
        return predict_passes(
            self.austin.lat, self.austin.lon, self.austin.h,
            self.iss_rv.rECEF, self.iss_rv.rECI, self.iss_rv.julian_date,
            min_elevation=self.min_elevation#, loc=location, sat=satellite
        )

    def mem_predict_passes(self):
        return predict_passes(
            self.austin.lat, self.austin.lon, self.austin.h,
            self.iss_rv.rECEF, self.iss_rv.rECI, self.iss_rv.julian_date,
            min_elevation=self.min_elevation#, loc=location, sat=satellite
        )


if __name__=="__main__":
    print('Create benchmark suite')
    suite = Suite()
    print('Set up benchmark')
    suite.setup()
    print('Run time_predict_passes')
    suite.time_predict_passes()
    print('Run mem_predict_passes')
    suite.mem_predict_passes()


