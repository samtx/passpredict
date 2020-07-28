# models.py
from dataclasses import dataclass
from functools import update_wrapper

import numpy as np

from .schemas import Location, Point, Overpass
from .constants import RAD2DEG
from .rotations.transform import ecef2sez
from .rotations.rotations import site_ECEF
from .timefn import jday2datetime


class reify(object):
    """From skyfield.descriptorlib"""
    def __init__(self, method):
        self.method = method
        update_wrapper(self, method)

    def __get__(self, instance, objtype=None):
        if instance is None:
            return self
        value = self.method(instance)
        instance.__dict__[self.__name__] = value
        return value


class Sun():
    def __init__(self):
        self.rECI = None
        self.rECEF = None
        self.time = None  # time object

    
# class Time():
#     __slots__ = ['jd', 'jd_utc1', 'tt', 'tt_utc1', 'datetime']
#     def __init__(self, jd=None, jd_utc1=None, tt=None, tt_utc1=None, datetime=None):
#         self.jd = jd
#         self.jd_utc1 = jd_utc1
#         self.tt = tt
#         self.tt_utc1 = tt_utc1
#         self.datetime = datetime       # datetimes

#     def __hash__(self):
#         return hash(str(t.jd))


class SpaceObject():
    # __slots__ = ['time', 'rECEF', 'rECI','latitude','longitude','altitude','illuminated', 'rECEF_astropy', 'subpoint']
    def __init__(self):
        self.time = None  # time object
        self.rECEF = None
        self.rECI = None
        self.latitude = None
        self.longitude = None
        self.altitude = None
        self.illuminated = None
        self.rECEF_astropy = None
        self.subpoint = None
        self.meta = None


class RhoVector():
    """
    Vector from topographic location to space object
    """
    # __slots__ = ['time', 'rSEZ', 'rECEF', 'rng', 'az', 'el', 'ra', 'dec', 'sat', 'location']
    def __init__(self, sat: SpaceObject, location: Location):
        self.sat = sat
        self.location = location
        self.time = sat.time

    @reify
    def rECEF(self):
        rsiteECEF = site_ECEF(self.location.lat, self.location.lon, self.location.h)
        return self.sat.rECEF - np.array([[rsiteECEF[0]],[rsiteECEF[1]],[rsiteECEF[2]]], dtype=np.float64)

    @reify
    def rSEZ(self):
        return ecef2sez(self.rECEF, self.location.lat, self.location.lon)

    @reify
    def rng(self):
        return np.linalg.norm(self.rSEZ, axis=0)

    @reify
    def el(self):
        return np.arcsin(self.rSEZ[2] / self.rng) * RAD2DEG

    @reify
    def az(self):
        tmp = np.arctan2(self.rSEZ[0], self.rSEZ[1])
        az = (tmp + np.pi * 0.5) * RAD2DEG
        idx = np.all([self.rSEZ[0] < 0, self.rSEZ[1] < 0], axis=0)
        az[idx] %= 360 
        return az

    def point(self, idx):
        return Point.construct(
            datetime=jday2datetime(self.time.jd[idx]),
            azimuth=self.az[idx],
            elevation=self.el[idx],
            range=self.rng[idx]
        )

    def find_overpasses(self, min_elevation=10, store_sat_id=False):
        # Find Overpasses
        el0 = self.el[:-1] - min_elevation
        el1 = self.el[1:] - min_elevation
        el_change_sign = (el0*el1 < 0)   
        start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]  # Find the start of an overpass
        end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]    # Find the end of an overpass
        num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
        if start_idx.size < end_idx.size:
            end_idx = end_idx[1:]
        sat_overpasses = [None] * num_overpasses
        for j in range(num_overpasses):
            # Store indecies of overpasses in a list
            idx0 = start_idx[j]
            idxf = end_idx[j]
            overpass_idx = np.arange(idx0, idxf+1, dtype=int)
            idxmax = np.argmax(self.el[overpass_idx])
            start_pt = self.point(idx0)
            max_pt = self.point(idx0 + idxmax)
            end_pt = self.point(idxf)
            if store_sat_id:
                overpass = Overpass.construct(
                    # satellite_id=self.meta.id,
                    start_pt=start_pt,
                    max_pt=max_pt,
                    end_pt=end_pt
                )
            else:
                overpass = Overpass.construct(
                    start_pt=start_pt,
                    max_pt=max_pt,
                    end_pt=end_pt
                )
            sat_overpasses[j] = overpass
        return sat_overpasses


class SatelliteRV():
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



