# models.py
from dataclasses import dataclass


import numpy as np


class Sun():
    def __init__(self):
        rECI = None
        rECEF = None
        t = None  # time object

    
class Time():
    def __init__(self):
        jd = None
        jd_utc = None
        tt = None
        tt_utc = None
        dt = None       # datetimes

    def __hash__(self):
        return hash(str(t.jd))


class SpaceObject():
    __slots__ = ['t', 'rECEF', 'rECI','latitude','longitude','altitude','visible']
    def __init__(self):
        self.t = None
        self.rECEF = None
        self.rECI = None
        self.latitude = None
        self.longitude = None
        self.altitude = None
        self.visible = None


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



