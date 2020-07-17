# models.py
from dataclasses import dataclass
from functools import update_wrapper

import numpy as np


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

    
class Time():
    __slots__ = ['jd', 'jd_utc1', 'tt', 'tt_utc1', 'datetime']
    def __init__(self, jd=None, jd_utc1=None, tt=None, tt_utc1=None, datetime=None):
        self.jd = jd
        self.jd_utc1 = jd_utc1
        self.tt = tt
        self.tt_utc1 = tt_utc1
        self.datetime = datetime       # datetimes

    def __hash__(self):
        return hash(str(t.jd))


class SpaceObject():
    __slots__ = ['time', 'rECEF', 'rECI','latitude','longitude','altitude','illuminated', 'rECEF_astropy', 'subpoint']
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


class RhoVector():
    """
    Vector from topographic location to space object
    """
    __slots__ = ['time', 'rSEZ', 'rECEF', 'rng', 'az','el','ra','dec']
    def __init__(self):
        self.time = None
        self.rSEZ = None
        self.rECEF = None
        self.rng = None
        self.az = None
        self.el = None
        self.ra = None
        self.dec = None
    

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



