# models
import datetime
import numpy as np


class Orbit():
    """Results from propagating a TLE"""
    def __init__(self, dt, tsince, rTEME, rECEF, tle_str, sat_id, dt_start):
         self.dt = dt
         self.tsince = tsince
         self.rTEME = rTEME
         self.rECEF = rECEF
         self.tle_str = tle_str
         self.sat_id = sat_id
         self.dt_start = dt_start
         self.modified = None

    def update_modified(self):
        self.modified = datetime.datetime.now()

class Point():
    """A point in the sky in topocentric horizon coordinates"""
    def __init__(self, dt, az, el, rng):
        self.dt = dt
        self.az = az
        self.el = el
        self.rng = rng

    def __repr__(self):
        dtstr = self.dt.strftime(r'%b %d %Y, %H:%M:%S')
        s = f'{dtstr}, az={self.az:.1f} deg, el={self.el:.1f} deg, rng={self.rng:.1f} km'
        return s


class Overpass():
    """An overpass of a satellite"""
    def __init__(self, start_point, max_point, end_point):
        self.start_point = start_point
        self.max_point = max_point
        self.end_point = end_point
        self.t = None
        self.r = None
        self.el = None
        self.az = None
        self.rng = None


def process_overpasses(overpasses, t, az, el, rng, dt_start):
    for j, overpass in enumerate(overpasses):
        overpass_len = len(overpass)
        start_idx = overpass[0]
        end_idx = overpass[overpass_len]
        start_pt = Point()


def tsince_to_datetime(tsince, dt_start):
    """Convert vector of minutes since epoch to an array of corresponding
    datetime values. Assumes that tsince is evenly spaced.
    Args:
        tsince : float (n)
        dt_start : datetime, beginning of epoch
    Returns:
        dt_ary : datetime (n), array of datetime objects
    """
    n = len(tsince)
    minute_step = datetime.timedelta(minutes=tsince[1] - tsince[0])
    dt_end = dt_start + datetime.timedelta(minutes=tsince[n-1])
    dt_ary = np.arange(dt_start, dt_end, minute_step)
    return dt_ary.astype(datetime.datetime)


def