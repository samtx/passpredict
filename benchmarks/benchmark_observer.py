# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import datetime

import numpy as np

from passpredict import *
from passpredict.locations import Location
from passpredict.satellites import SatellitePredictor
from passpredict import MemoryTLESource
from passpredict.observers import Observer, BruteForceObserver
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
        source = MemoryTLESource()
        source.add_tle(TLE(satid, tle_lines))
        satellite = SatellitePredictor(satid, source)

        location = Location("Austin, Texas", 30.2711, -97.7434, 0)
        self.min_elevation = 10
        self.observer = passpredict.Observer(location, satellite, aos_at_dg=self.min_elevation)
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

    def track_elevation_at_function_calls(self):
        passpredict.predict_single_satellite_overpasses(
            self.satellite,
            self.location,
            self.start,
            self.ndays,
            10,  # min elevation
        )
        res = self.observer._elevation_at.cache_info()
        return res.hits + res.misses

    def track_elevation_at_function_cache_ration(self):
        passpredict.predict_single_satellite_overpasses(
            self.satellite,
            self.location,
            self.start,
            self.ndays,
            10,  # min elevation
        )
        res = self.observer._elevation_at.cache_info()
        return res.hits / (res.misses + res.hits)



class PredictOverpassesBruteForce:
    """
    ISS prediction over Austin, Texas with Brute Force algorithm
    """
    params = (
        [60, 20, 5],
        [0.5, 0.25, 0.1],
    )
    param_names = ['time_step', 'tolerance_s']

    def setup(self, time_step, tolerance_s):
        self.start = datetime.datetime(2020, 7, 14, tzinfo=datetime.timezone.utc)
        self.ndays = 10
        self.end = self.start + datetime.timedelta(days=self.ndays)
        satid = 'ISS'
        tle_lines = (
            "1 25544U 98067A   20196.51422950 -.00000046  00000-0  72206-5 0  9999",
            "2 25544  51.6443 213.2207 0001423 114.8006 342.8278 15.49514729236251"
        )
        source = MemoryTLESource()
        source.add_tle(TLE(satid, tle_lines))
        satellite = SatellitePredictor(satid, source)

        location = Location("Austin, Texas", 30.2711, -97.7434, 0)
        self.min_elevation = 10
        self.observer = passpredict.BruteForceObserver(location, satellite, aos_at_dg=self.min_elevation)

    def time_brute_force_observer(self, *args):
        pass_iterator = self.observer.iter_passes(self.start, limit_date=self.end)
        list(pass_iterator)