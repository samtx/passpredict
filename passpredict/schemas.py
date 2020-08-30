# schemas
import datetime
from typing import List
from enum import Enum
from functools import cached_property

import numpy as np
from pydantic import BaseModel, Field

from .timefn import jday2datetime
from .tle import TleSchema as Tle

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
    
    @classmethod
    def from_rho(cls, rho, idx):
        """Create a Point object directly from the rho vector and index without validation"""
        return cls.construct(
            datetime=jday2datetime(rho.time.jd[idx]),
            azimuth=rho.az[idx],
            elevation=rho.el[idx],
            range=rho.rng[idx]
        )


    # def __repr__(self):
    #     dtstr = self.datetime.strftime("%b %d %Y, %H:%M:%S")
    #     s = "{}UTC el={:.1f}d, az={:.1f}d, rng={:.1f}km".format(
    #         dtstr, self.elevation, self.azimuth, self.range)
    #     return s


class Location(BaseModel):
    # __slots__ = ['lat', 'lon', 'h', 'name', 'tz']
    lat: float       # latitude, decimal degrees, positive is North
    lon: float       # longitude, decimal degrees, positive is East
    h: float = 0.0   # elevation [m]
    name: str = None
    # tz: Timezone = None  # timezone object


class Satellite(BaseModel):
    id: int    # NORAD ID
    name: str = None


class SatelliteDetails(BaseModel):
    id = int   # NORAD ID
    name: str = None
    cospar_id: str = None
    radio_downlink: float = None  # MHz
    intrinsic_brightness: float = None  # magnitude
    maximum_brightness: float = None    # magnitude
    category: str = None
    description: str = None
    country: str = None
    launch_date: datetime.date = None
    launch_site: str = None
    mass: float = None   # kg


class PassType(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'


class Overpass(BaseModel):
    start_pt: Point
    max_pt: Point
    end_pt: Point
    satellite_id: int = None
    type: PassType = None
    vis_start_pt: Point = None
    vis_end_pt: Point = None
    brightness: float = None
    

class OverpassResultBase(BaseModel):
    location: Location


class SingleSatOverpassResult(OverpassResultBase):
    satellite: Satellite
    overpasses: List[Overpass]


class MultiSatOverpassResult(OverpassResultBase):
    overpasses: List[Overpass]

