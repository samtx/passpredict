# schemas
import datetime
from dataclasses import dataclass
from typing import Any, List
from enum import Enum

import numpy as np
from pydantic import BaseModel

# From Pydantic, to use Numpy arrays
# ref: https://github.com/samuelcolvin/pydantic/issues/380#issuecomment-620378743
class _ArrayMeta(type):
    def __getitem__(self, t):
        return type('Array', (Array,), {'__dtype__': t})

class Array(np.ndarray, metaclass=_ArrayMeta):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate_type

    @classmethod
    def validate_type(cls, val):
        dtype = getattr(cls, '__dtype__', None)
        if isinstance(dtype, tuple):
            dtype, shape = dtype
        else:
            shape = tuple()

        result = np.array(val, dtype=dtype, copy=False, ndmin=len(shape))
        assert not shape or len(shape) == len(result.shape)  # ndmin guarantees this

        if any((shape[i] != -1 and shape[i] != result.shape[i]) for i in range(len(shape))):
            result = result.reshape(shape)
        return result

"""
Example use: 

class Model(pydantic.BaseModel):
    int_values: Array[float]
    any_values: Array
    shaped1_values: Array[float, (-1, )]
    shaped2_values: Array[float, (2, 1)]
    shaped3_values: Array[float, (4, -1)]
    shaped4_values: Array[float, (-1, 4)]
"""

# Create timezone field
# ref: https://pydantic-docs.helpmanual.io/usage/types/#enums-and-choices
# ref: https://pydantic-docs.helpmanual.io/usage/types/#custom-data-types
# class Timezone(pytz.timezone):
#     @classmethod
#     def __get_validators__(cls):
#         yield cls.validate

#     @classmethod
#     def validate(cls, v):
#         if not isinstance(v, pytz.timezone):
#             raise TypeError('not a pytz timezone')
        
#     def __repr__(self):
#         return super().__repr__()


class Timezone(BaseModel):
    offset: float  # UTC offset
    name: str = None
    def __repr__(self):
        return name


class Sun(BaseModel):
    rECI: Array[float]
    jdt: Array[float]
    dt_ary: Array[float]


COORDINATES = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','N']


class Point(BaseModel):
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

# utc = Timezone(offset=0.0, name='UTC')

class Location(BaseModel):
    # __slots__ = ['lat', 'lon', 'h', 'name', 'tz']
    lat: float       # latitude, decimal degrees, positive is North
    lon: float       # longitude, decimal degrees, positive is East
    h: float = 0.0   # elevation [m]
    name: str = None
    # tz: Timezone = None  # timezone object


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


class Satellite(BaseModel):
    id: int
    name: str = None


class Tle(BaseModel):
    # __slots__ = ['tle1','tle2','epoch','satellite']
    tle1: str
    tle2: str
    epoch: datetime.datetime
    satellite: Satellite


class Overpass(BaseModel):
    start_pt: Point
    max_pt: Point
    end_pt: Point
    satellite_id: int = None
    

class OverpassResultBase(BaseModel):
    location: Location


class SingleSatOverpassResult(OverpassResultBase):
    satellite: Satellite
    overpasses: List[Overpass]


class MultiSatOverpassResult(OverpassResultBase):
    overpasses: List[Overpass]


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