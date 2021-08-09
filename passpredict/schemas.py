# schemas
import datetime
from typing import List
from enum import Enum
from functools import cached_property

import numpy as np
from pydantic import BaseModel, Field

from .timefn import jday2datetime, epoch_to_jd
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


class Orbit(BaseModel):
    satnum: str
    jdsatepoch: float
    jdsatepochF: float
    bstar: float
    ecco: float
    inclo: float
    nodeo: float
    argpo: float
    mo: float
    no_kozai: float
    revnum: int
    elnum: int
    classification: str
    ephtype: int

    @classmethod
    def from_tle(cls, tle1, tle2):
        """
        Convert TLE strings to Orbit object
        """
        satnum = tle1[2:7]
        classification = tle1[8]
        epoch_year = int(tle1[18:20])
        epoch_days = float(tle1[20:32])
        jdsatepoch, jdsatepochF = epoch_to_jd(epoch_year, epoch_days)
        # ndot = float(tle1[34:44])
        # nddot = float(tle1[45:52])
        bstar = float(tle1[54:62])
        ephtype = tle1[63]
        elnum = int(tle1[65:69])
        inclo = float(tle2[9:17])  # inclination
        nodeo = float(tle2[18:26])  # right ascension of ascending node
        ecco = float(tle2[27:34]) / 1e7  # eccentricity
        argpo = float(tle2[35:43])
        mo = float(tle2[44:52])    # mean anomaly
        no_kozai = float(tle2[53:64])   # mean motion
        revnum = int(tle2[64:69])
        orbit = cls(
            satnum=satnum,
            jdsatepoch=jdsatepoch,
            jdsatepochF=jdsatepochF,
            bstar=bstar,
            inclo=inclo,
            nodeo=nodeo,
            ecco=ecco,
            argpo=argpo,
            mo=mo,
            no_kozai=no_kozai,
            revnum=revnum,
            elnum=elnum,
            classification=classification,
            ephtype=ephtype
        )
        return orbit

    @classmethod
    def from_omm(cls, omm):
        """
        Convert OMM object to Orbit object
        """
        orbit = cls(
            satnum=omm.satnum,
            jdsatepoch=omm.jdsatepoch,
            jdsatepochF=omm.jdsatepochF,
            bstar=omm.bstar,
            inclo=omm.inclo,
            nodeo=omm.nodeo,
            ecco=omm.ecco,
            argpo=omm.argpo,
            mo=omm.mo,
            no_kozai=omm.no_kozai,
            revnum=omm.revnum,
            elnum=omm.elnum,
            classification=omm.classification,
            ephtype=omm.ephtype
        )
        return orbit



class Satellite(BaseModel):
    id: int    # NORAD ID
    name: str = None


class SatelliteCategory(str, Enum):
    military = 'military'
    communication = 'communication'
    earth_observing = 'earth observing'
    space_science = 'space science'



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
    perigee: float = None  # km
    apogee: float = None   # km


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

