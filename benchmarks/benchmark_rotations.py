# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.

from math import radians
import numpy as np

from passpredict import _rotations


class Rotations:
    """
    Example from Vallado, Eg 11-6, p. 912
    """
    def setup(self, *args):
        self.lat = radians(42.38)
        self.lon = radians(-71.13)
        self.location_ecef = np.array([1526.122, -4465.064, 4276.894])
        self.satellite_ecef = np.array([885.7296, -4389.3856, 5070.1765])

    def time_ecef_to_razel(self):
        _rotations.razel(self.lat, self.lon, self.location_ecef, self.satellite_ecef)

    def time_elevation_at(self):
        _rotations.elevation_at(self.lat, self.lon, self.location_ecef, self.satellite_ecef)

    def time_range_at(self):
        _rotations.range_at(self.lat, self.lon, self.location_ecef, self.satellite_ecef)

    def time_ecef_to_llh(self):
        _rotations.ecef_to_llh(self.satellite_ecef)


class SolarRotations:
    """
    Example from Vallado, Eg.5-1, p.280, April 2, 2006, 00:00 UTC
    """

    def setup(self, *args):
        self.jd = 2453827.5
        self.rmod = np.array([146186212.0, 28788976.0, 12481064.0])
        self.rpef = np.empty(3, dtype=np.double)

    def time_mod2ecef(self):
        _rotations.mod2ecef(self.jd, self.rmod, self.rpef)
