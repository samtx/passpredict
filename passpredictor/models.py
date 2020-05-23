# models
import datetime
import numpy as np
from dataclasses import dataclass
import pytz

class Sun(object):
    def __init__(self):
        self.rECI = None
        self.jdt = None
        self.dt_ary = None

COORDINATES = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','N']

@dataclass
class Point:
    # __slots__ = ['datetime', 'azimuth', 'elevation', 'range', 'declination', 'right_ascension']
    datetime: datetime.datetime
    azimuth: float
    elevation: float
    range: float
    declination: float = None
    right_ascension: float = None

    def direction_from_azimuth(self):
        ''' Return direction from azimuth degree '''
        azm = self.azimuth % 360
        mod = 360/16. # number of degrees per coordinate heading
        start = 0 - mod/2
        n = np.floor((azm-start)/mod).astype(int)
        return COORDINATES[n]


    # def __repr__(self):
    #     dtstr = self.datetime.strftime("%b %d %Y, %H:%M:%S")
    #     s = "{}UTC el={:.1f}d, az={:.1f}d, rng={:.1f}km".format(
    #         dtstr, self.elevation, self.azimuth, self.range)
    #     return s


class Overpass(object):
    __slots__ = ['location', 'satellite', 'start_pt', 'max_pt', 'end_pt', 't', 'r']
    def __init__(self, location, satellite, start_pt, max_pt, end_pt, t, r):
        self.location = location
        self.satellite = satellite
        self.start_pt = start_pt
        self.max_pt = max_pt
        self.end_pt = end_pt
        self.t = t
        self.r = r



@dataclass
class Location:
    # __slots__ = ['lat', 'lon', 'h', 'name', 'tz']
    lat: float       # latitude, decimal degrees, positive is North
    lon: float       # longitude, decimal degrees, positive is East
    h: float = 0.0   # elevation [m]
    name: str = None
    tz: pytz.tzinfo = pytz.utc  # timezone object


class SatelliteRV(object):
    __slots__ = ['satellite','tle','rsun','datetime','julian_date','rECEF',
                 'rECI','latitude','longitude','altitude','visible']
    def __init__(self):
        self.satellite = None
        self.tle = None
        self.rsun = None
        self.datetime = None
        self.julian_date = None
        self.rECEF = None
        self.rECI = None
        self.latitude = None
        self.longitude = None
        self.altitude = None
        self.visible = None

class Tle(object):
    __slots__ = ['tle1','tle2','epoch','satellite']
    def __init__(self, tle1, tle2, epoch, satellite):
        self.tle1 = tle1
        self.tle2 = tle2
        self.epoch = epoch
        self.satellite = satellite


@dataclass
class Satellite:
    id: int
    name: str


def process_overpasses(overpasses, t, az, el, rng, dt_start):
    for j, overpass in enumerate(overpasses):
        overpass_len = len(overpass)
        start_idx = overpass[0]
        end_idx = overpass[overpass_len]
        start_pt = Point()


def tsince_to_datetime(tsince, dt_start):
    """
    Convert vector of minutes since epoch to an array of corresponding
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