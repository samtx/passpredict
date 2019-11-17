# models
import datetime
import numpy as np


class Point(object):
    def __init__(self, datetime, azimuth, elevation, range_):
        self.datetime = datetime
        self.azimuth = azimuth
        self.elevation = elevation
        self.range = range_
    def __repr__(self):
        dtstr = self.datetime.strftime("%b %d %Y, %H:%M:%S")
        s = "{}UTC el={:.1f}d, az={:.1f}d, rng={:.1f}km".format(
            dtstr, self.elevation, self.azimuth, self.range)
        return s


class Overpass(object):
    def __init__(self, start_pt, max_pt, end_pt, t, r):
        self.start_pt = start_pt
        self.max_pt = max_pt
        self.end_pt = end_pt
        self.t = t
        self.r = r
        self.sat = None
        self.location = None


class Location(object):
    def __init__(self, lat, lon, h, name=None, tz=None):
        self.lat = lat
        self.lon = lon
        self.h = h
        self.name = name
        self.tz = tz
    def __repr__(self):
        return self.name


class SatelliteRV(object):
    def __init__(self):
        self.satellite = None
        self.tle = None
        self.rsun = None
        self.dt = None
        self.jdt = None
        self.rECEF = None
        self.rECI = None
        self.modified = None
        self.deltaT = None
        self.lat = None
        self.lon = None
        self.alt = None
        self.is_illum = None


class SunPosition(object):
    def __init__(self):
        self.rECI = None
        self.jdt = None
        self.dt_ary = None


class Tle(object):
    def __init__(self, tle1, tle2, dt):
        self.tle1 = tle1
        self.tle2 = tle2
        self.dt = dt
        self.satellite = None


class SatelliteDetail(object):
    def __init__(self, satid, name):
        self.id = satid
        self.name = name


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