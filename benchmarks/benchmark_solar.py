# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.

import numpy as np

from passpredict import solar
from passpredict import _solar


class SunPosition:
    def setup(self):
        self.jd = 2450540.547222222
        self.rmod = np.empty(3, dtype=np.double)

    def time_sun_pos_cython(self):
        _solar.sun_pos(self.jd)

    def time_sun_pos(self):
        solar.sun_pos.__wrapped__(self.jd)

    def time_sun_pos_mod(self):
        _solar.sun_pos_mod(self.jd, self.rmod)


class SunPositionCache:
    def setup(self):
        self.jd = 2450540.547222222
        self.rmod = np.empty(3, dtype=np.double)
        solar.sun_pos(self.jd)

    def time_sun_pos_cache(self):
        solar.sun_pos(self.jd)


class SatelliteIllumination:
    def setup(self):
        self.rsat = np.array([885.7296, -4389.3856, 5070.1765])
        self.rsun = np.array([-1.43169570e+08, 4.12567046e+07, 1.27096677e+07])

    def time_sat_illumination_distance(self):
        _solar.sat_illumination_distance(self.rsat, self.rsun)

