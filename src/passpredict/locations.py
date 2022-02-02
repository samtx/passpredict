from functools import cached_property
from datetime import datetime, timezone as py_timezone
from math import radians

import numpy as np
from orbit_predictor.locations import Location as LocationBase
from orbit_predictor import coordinate_systems

from .zoneinfo import ZoneInfo
from .utils import get_timezone_from_latlon
from .time import julian_date_from_datetime
from .solar import sun_pos
from ._rotations import elevation_at


class Location(LocationBase):

    def __init__(self, name, latitude_deg, longitude_deg, elevation_m):
        """Location.
        Parameters
        ----------
        latitude_deg : float
            Latitude in degrees.
        longitude_deg : float
            Longitude in degrees.
        elevation_m : float
            Elevation in meters.
        """
        self.name = name
        self.latitude_deg = latitude_deg
        self.longitude_deg = longitude_deg
        self.latitude_rad = radians(latitude_deg)
        self.longitude_rad = radians(longitude_deg)
        self.elevation_m = elevation_m
        self.position_ecef = coordinate_systems.geodetic_to_ecef(
            self.latitude_rad,
            self.longitude_rad,
            elevation_m / 1000.)
        self.position_llh = latitude_deg, longitude_deg, elevation_m
        self.recef = np.array(self.position_ecef)

    def dict(self) -> dict:
        d = {
            'name': self.name,
            'lat': self.lat,
            'lon': self.lon,
            'h': self.h
        }
        return d

    @property
    def lat(self) -> float:
        return self.latitude_deg

    @property
    def lon(self) -> float:
        return self.longitude_deg

    @property
    def h(self) -> float:
        return self.elevation_m

    @cached_property
    def timezone(self) -> ZoneInfo:
        """ Find timezone """
        return get_timezone_from_latlon(self.latitude_deg, self.longitude_deg)

    @property
    def tz(self) -> ZoneInfo:
        return self.timezone

    @cached_property
    def offset(self) -> float:
        """  Compute timezone offset in hours from UTC  """
        now = datetime.now(self.timezone)
        delta = now.utcoffset().total_seconds() / 3600
        return delta

    def sun_elevation_jd(self, jd: float) -> float:
        """
        Computes elevation angle of sun relative to location. Returns degrees.
        """
        sun_recef = sun_pos(jd)
        el = elevation_at(self.latitude_rad, self.longitude_rad, self.recef, sun_recef)
        return el

    def sun_elevation(self, dt: datetime) -> float:
        """
        Computes elevation angle of sun relative to location. Returns degrees.
        """
        jd, jdfr = julian_date_from_datetime(dt)
        jd = jd + jdfr
        el = self.sun_elevation_jd(jd)
        return el

    def is_sunlit(self, dt: datetime) -> bool:
        """
        Computes elevation angle of sun relative to location
        Returns True if elevation > -6 degrees
        """
        el = self.sun_elevation(dt)
        return el > -6

    def __repr__(self):
        deg = u'\N{DEGREE SIGN}'
        s = '<Location '
        if self.name:
            s += self.name + ' '
        s += f'({self.latitude_deg}{deg} , {self.longitude_deg}{deg})'
        s += '>'
        return s


