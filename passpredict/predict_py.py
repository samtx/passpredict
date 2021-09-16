# Python only module for prediction functions
from functools import partial
import datetime

from passpredict import Location, Satellite
from passpredict.predict import compute_elevation_angle
from passpredict.timefn import julian_date, jday2datetime


class Observer:
    def __init__(self, location: Location, satellite: Satellite):
        self.location = location
        self.satellite = satellite

    def compute_elevation(self, jd: float) -> float:
        """
        Compute elevation angle at julian date jd
        """
        el = compute_elevation_angle(jd, self.location, self.satellite)
        return el

    def predict_overpasses(self, jd0: float = None, jdmax: float = None):
        if jd0 is None:
            jd0 = self.satellite.epoch
        if jdmax is None:
            jdmax = jd0 + 10

        # get overpasses over time frame
        overpasses = []
        jd = jd0
        while (jd < jdmax):
            jd_aos = self.find_next_aos(jd, jdmax)
            overpasses.append(jday2datetime(jd_aos))
            if jd_aos < 0:
                break
            # find next LOS
            jd = jd_aos + 5 / 86400.0  # add 5 seconds
            jd_los = self.find_next_los(jd, jdmax)
            if jd_los < 0:
                break
            jd = jd_los + 20 / 1440.0  # advance 20 minutes
        return overpasses

    def find_next_aos(self, jd0, jdf):
        el = -90
        # coarse time steps
        dt_large = 30 / 86400.0  # 30 seconds
        jd = jd0 - dt_large
        k = 0
        while (el < 0.0) and (jd <= jdf):
            jd += dt_large
            el = self.compute_elevation(jd)
            # k += 1
        if jd > jdf:
            return -1
        # fine time steps
        dt_small = 1 / 86400.0  # 1 second
        jd = jd - dt_large
        j = 0
        el = self.compute_elevation(jd)
        while (el < 0) and (jd <= jdf):
            jd += dt_small
            el = self.compute_elevation(jd)
            # j += 1
        if jd > jdf:
            return -1
        # print(f'k={k}  j={j}')
        return jd

    def find_next_los(self, jd0, jdf):
        """
        Find next loss of signal from time jd0
        """
        el = +90
        # coarse time steps
        dt_large = 20 / 86400.0  # 20 seconds
        jd = jd0 - dt_large
        k = 0
        while (el > 0.0) and (jd <= jdf):
            jd += dt_large
            el = self.compute_elevation(jd)
            # k += 1
        if jd > jdf:
            return -1
        # fine time steps
        dt_small = 1 / 86400.0  # 1 second
        jd = jd - dt_large
        j = 0
        el = self.compute_elevation(jd)
        while (el > 0) and (jd <= jdf):
            jd += dt_small
            el = self.compute_elevation(jd)
            # j += 1
        if jd > jdf:
            return -1
        # print(f'k={k}  j={j}')
        return jd


def predict(jd0, jdf, location, satellite):
    """
    Main prediction algorithm for finding overpasses
    """
    observer = Observer(location=location, satellite=satellite)
    res = observer.predict_overpasses()
    return [jday2datetime(jd)]


