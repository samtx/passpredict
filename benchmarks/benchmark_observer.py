# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import datetime

import numpy as np

from passpredict import *
from passpredict.locations import Location
from passpredict.satellites import SatellitePredictor
from passpredict import MemoryTLESource, Observer
import passpredict



class PredictOverpasses:
    """
    ISS prediction over Austin, Texas
    """
    def setup(self, *args):
        self.start = datetime.datetime(2020, 7, 14, tzinfo=datetime.timezone.utc)
        self.ndays = 10
        self.end = self.start + datetime.timedelta(days=self.ndays)
        satid = 'ISS'
        tle_lines = (
            "1 25544U 98067A   20196.51422950 -.00000046  00000-0  72206-5 0  9999",
            "2 25544  51.6443 213.2207 0001423 114.8006 342.8278 15.49514729236251"
        )
        satellite = SatellitePredictor(satid)
        satellite.tle = TLE(satid, tle_lines, self.start)
        satellite.set_propagator()

        location = Location("Austin, Texas", 30.2711, -97.7434, 0)
        self.observer = passpredict.Observer(location, satellite)
        self.satellite = satellite
        self.location = location

    def time_observer_iter_passes(self):
        pass_iterator = self.observer.iter_passes(self.start, limit_date=self.end)
        list(pass_iterator)

    def time_predict_single_satellite_overpass(self):
        passpredict.predict_single_satellite_overpasses(
            self.satellite,
            self.location,
            self.start,
            self.ndays,
            10,  # min elevation
        )
