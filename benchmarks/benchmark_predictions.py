# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import datetime

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from passpredict.predictions import predict, find_overpasses, compute_sun_data, compute_time_array, compute_satellite_data
from passpredict.schemas import Location, Satellite, Tle, Point, Overpass
from passpredict.models import SpaceObject, RhoVector, Sun, Sat
from passpredict.propagate import propagate_satellite
from passpredict.timefn import julian_date, jd2jc, jd2utc1, jday2datetime
from passpredict.rotations.polar import eop
from passpredict.rotations.rotations import site_ECEF
from passpredict.rotations.transform import ecef2eci, ecef2eci, ecef2sez, teme2ecef, teme2eci
from passpredict.solar import sun_pos, is_sat_illuminated
from passpredict.topocentric import razel, site_sat_rotations
from passpredict.utils import epoch_from_tle

from .utils import save_pickle, load_pickle


class Predict:
    """
    An example benchmark that times the performance of various kinds
    of iterating over dictionaries in Python.
    """

    # SECONDS_PER_DAY = 3600*24*7

    # self.params = ([SECONDS_PER_DAY*i for i in [7, 10, 14, 18, 21]])
    # self.param_names = ['num_time_steps']

    def setup(self, *args):
        satellite = Satellite(id=25544, name='ISS')
        location = Location(lat=30.2711, lon=-97.7434, h=0, name='Austin, Texas')
        dt_seconds = 1.0
        tle1 = '1 25544U 98067A   20196.51422950 -.00000046  00000-0  72206-5 0  9999'
        tle2 = '2 25544  51.6443 213.2207 0001423 114.8006 342.8278 15.49514729236251'    
        tle = Tle(
            tle1=tle1,
            tle2=tle2,
            epoch=epoch_from_tle(tle1),
            satellite=satellite
        ) 
        dt_start = datetime.datetime(2020, 7, 14, 11, 17, 00, tzinfo=datetime.timezone.utc)
        dt_end = dt_start + datetime.timedelta(days=14)
        t = compute_time_array(dt_start, dt_end, dt_seconds)
        sun = compute_sun_data(t)
        sat = compute_satellite_data(tle, t, sun)
        self.location = location
        self.sat = sat
        self.t = t
        self.sun = sun
        
        
    def time_find_overpasses(self):
        out = find_overpasses(self.location, [self.sat], self.t, self.sun)


    def peakmem_predict_realtime_compute(self):
        out = find_overpasses(self.location, [self.sat], self.t, self.sun)