# models.py
from dataclasses import dataclass
from functools import update_wrapper

import numpy as np

from .schemas import Location, Point, Overpass
from .constants import RAD2DEG
from .rotations import ecef2sez
from .topocentric import site_ECEF
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


class SpaceObject:
    # __slots__ = ['time', 'rECEF', 'rECI','latitude','longitude','altitude','illuminated', 'rECEF_astropy', 'subpoint']
    def __init__(self):
        self.time = None  # time object
        self.rECEF = None
        self.rECI = None
        
        
class Sun(SpaceObject):
    pass


class Sat(SpaceObject):
    def __init__(self):
        self.latitude = None
        self.longitude = None
        self.altitude = None
        self.illuminated = None
        self.rECEF_astropy = None
        self.subpoint = None
        self.id = None


class RhoVectorBase:
    pass


class RhoVector():
    """
    Vector from topographic location to space object
    """
    # __slots__ = ['time', 'rSEZ', 'rECEF', 'rng', 'az', 'el', 'ra', 'dec', 'sat', 'location']
    def __init__(self, sat: SpaceObject, location: Location, sun: Sun = None):
        self.sat = sat
        self.location = location
        self.time = sat.time
        
        if sun is not None:
            # If the sun variable is set, it's time object must be identical to the sat
            # assert np.all(sun.time.jd[[0, -1]], sat.time.jd[[0, -1]])
            assert sun.time.jd[0] == sat.time.jd[0]
            assert sun.time.jd[-1] == sat.time.jd[-1]
            assert sat.illuminated is not None
            self.sun = sun
            self.site_sun_rho = RhoVector(sun, location)
        else:
            self.sun = sun
            self.site_sun_rho = None


    @reify
    def rsiteECEF(self):
        return site_ECEF(self.location.lat, self.location.lon, self.location.h)
    
    @reify
    def rECEF(self):
        return self.sat.rECEF - np.array([[self.rsiteECEF[0]],[self.rsiteECEF[1]],[self.rsiteECEF[2]]], dtype=np.float64)

    def _ecef2sez(self, idx):
        pass 

    @reify
    def rSEZ(self):
        return ecef2sez(self.rECEF, self.location.lat, self.location.lon)

    @reify
    def rng(self):
        return np.linalg.norm(self.rSEZ, axis=0)

    def _elevation(self, rZ, rng):
        return np.arcsin(self.rSEZ[2] / self.rng) * RAD2DEG

    @reify
    def el(self):
        return self._elevation(self.rSEZ[2], self.rng)

    def az(self, idx):
        rS = self.rSEZ[0, idx]
        rE = self.rSEZ[1, idx]
        tmp = np.arctan2(rS, rE)
        az = (tmp + np.pi * 0.5) * RAD2DEG
        if rS < 0 and rE < 0:
            az %= 360 
        # idx = np.all([self.rSEZ[0] < 0, self.rSEZ[1] < 0], axis=0)
        # az[idx] %= 360 
        return az

    def point(self, idx):
        return Point.construct(
            datetime=jday2datetime(self.time.jd[idx]),
            azimuth=self.az(idx),
            elevation=self.el[idx],
            range=self.rng[idx]
        )

    def _start_end_index(self, x):
        """
        Finds the start and end indecies when a 1D array crosses zero
        """
        x0 = x[:-1]
        x1 = x[1:]
        x_change_sign = (x0*x1 < 0)   
        start_idx = np.nonzero(x_change_sign & (x0 < x1))[0]  # Find the start of an overpass
        end_idx = np.nonzero(x_change_sign & (x0 > x1))[0]    # Find the end of an overpass
        return start_idx, end_idx
        

    def find_overpasses(self, min_elevation=10, store_sat_id=False, sunset_el=-6):
        start_idx, end_idx = self._start_end_index(self.el - min_elevation)
        num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
        if start_idx.size < end_idx.size:
            end_idx = end_idx[1:]
        sat_overpasses = [None] * num_overpasses
        for j in range(num_overpasses):
            # Store indecies of overpasses in a list
            idx0 = start_idx[j]
            idxf = end_idx[j]
            idxmax = np.argmax(self.el[idx0:idxf+1])
            start_pt = self.point(idx0)
            max_pt = self.point(idx0 + idxmax)
            end_pt = self.point(idxf)

            # Find visible start and end times
            if self.sun is not None:
                sun_sez = ecef2sez(self.sun.rECEF[:,idx0:idxf+1], self.location.lat, self.location.lon)
                sun_rng = np.linalg.norm(sun_sez, axis=0)
                sun_el = np.arcsin(sun_sez[2] / sun_rng) * RAD2DEG
                site_in_sunset = sun_el - sunset_el < 0
                site_in_sunset_idx = np.nonzero(site_in_sunset)[0]
                if site_in_sunset_idx.size == 0:
                    # site is always sunlit, so overpass is in daylight
                    visibility = 1
                else:
                    # get satellite illumination values for this overpass
                    sat_visible = (self.sat.illuminated[idx0:idxf+1] * site_in_sunset)
                    if np.any(sat_visible):
                        visibility = 3 # site in night, sat is illuminated
                        sat_visible_idx = np.nonzero(sat_visible)[0]
                        sat_visible_start_idx = sat_visible_idx.min()
                        sat_visible_end_idx = sat_visible_idx.max()
                        vis_start_pt = self.point(idx0 + sat_visible_start_idx)
                        vis_end_pt = self.point(idx0 + sat_visible_end_idx)
                    else:
                        visibility = 2 # nighttime, not illuminated (radio night)
            else:
                visibility = None
            overpass_dict = {
                'start_pt': start_pt,
                'max_pt': max_pt,
                'end_pt': end_pt,
            }
            if store_sat_id:
                overpass_dict['satellite_id'] = self.sat.id
            if (visibility is not None) and (visibility >= 3):
                overpass_dict['vis_start_pt'] = vis_start_pt
                overpass_dict['vis_end_pt'] = vis_end_pt
            overpass_dict['visibility'] = visibility
            overpass = Overpass.construct(**overpass_dict)
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



