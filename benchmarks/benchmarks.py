# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import datetime

import numpy as np

from passpredict.predictions import predict
from passpredict.schemas import Location, Satellite, Tle, Point, Overpass
from passpredict.models import SpaceObject, RhoVector, Time, Sun
from passpredict.propagate import propagate_satellite
from passpredict.timefn import julian_date, jd2jc, jd2utc1, jday2datetime
from passpredict.rotations.polar import eop
from passpredict.rotations.rotations import site_ECEF
from passpredict.rotations.transform import ecef2eci, ecef2eci, ecef2sez, teme2ecef, teme2eci
from passpredict.solar import sun_pos, is_sat_illuminated
from passpredict.topocentric import razel, site_sat_rotations
from passpredict.utils import epoch_from_tle


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

        jdt0 = julian_date(dt_start)
        jdtf = julian_date(dt_end)
        total_days = (dt_start-dt_end).total_seconds()/60
        dt_days = 1.0/(24*60*60.0)  # dt_seconds = 1
        t = Time(jd=np.arange(jdt0, jdtf, dt_days, dtype=float))
        t.tt = jd2jc(t.jd)
        dUTC1, xp, yp = eop(t.jd)
        t.jd_utc1 = t.jd + dUTC1
        t.tt_utc1 = jd2jc(t.jd_utc1)

        sat = SpaceObject()
        sat.time = t
        rTEME, _ = propagate_satellite(tle.tle1, tle.tle2, t.jd)
        sat.rECEF = teme2ecef(rTEME, t.jd_utc1, xp, yp)
        sat.rECI = rTEME.view()

        # Compute sun-satellite quantities
        sun = Sun()
        sun.time = t
        sun.rECI = sun_pos(sun.time.jd)  # to do: use cached value
        sat.illuminated = is_sat_illuminated(rTEME.copy(), sun.rECI)


        self.location = location
        self.t = t
        self.sat = sat
        self.sun = sun

    def time_predict_realtime_compute(self):
        location = self.location
        t = self.t
        sat = self.sat
        sun = self.sun

        store_sat_id = False
        
        min_elevation = 10.1

        rho = RhoVector()
        rho.time = t
        rsiteECEF = site_ECEF(location.lat, location.lon, location.h)
        rho.rECEF = site_sat_rotations(rsiteECEF, sat.rECEF)
        rho.rSEZ = ecef2sez(rho.rECEF, location.lat, location.lon)
        rng, az, el = razel(rho.rSEZ)
        rho.rng = rng
        rho.az = az
        rho.el = el
        # rsiteECI = ecef2eci(rsiteECEF, jdt_utc1)
        
        # Find Overpasses
        el0 = rho.el[:-1] - min_elevation
        el1 = rho.el[1:] - min_elevation
        el_change_sign = (el0*el1 < 0)   
        start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]  # Find the start of an overpass
        end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]    # Find the end of an overpass
        num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
        if start_idx.size < end_idx.size:
            end_idx = end_idx[1:]
        overpasses = [None] * num_overpasses
        for j in range(num_overpasses):
            # Store indecies of overpasses in a list
            idx0 = start_idx[j]
            idxf = end_idx[j]
            overpass_idx = np.arange(idx0, idxf+1, dtype=int)
            idxmax = np.argmax(el[overpass_idx])
            start_pt = Point(
                datetime=jday2datetime(rho.time.jd[idx0]),
                azimuth=rho.az[idx0],
                elevation=rho.el[idx0],
                range=rho.rng[idx0]
            )
            max_pt = Point(
                datetime=jday2datetime(rho.time.jd[idx0 + idxmax]),
                azimuth=rho.az[idx0 + idxmax],
                elevation=rho.el[idx0 + idxmax],
                range=rho.rng[idx0 + idxmax]
            )
            end_pt = Point(
                datetime=jday2datetime(rho.time.jd[idxf]),
                azimuth=rho.az[idxf],
                elevation=rho.el[idxf],
                range=rho.rng[idxf]
            )
            if store_sat_id:
                overpass = Overpass(
                    satellite_id=satellite.id,
                    start_pt=start_pt,
                    max_pt=max_pt,
                    end_pt=end_pt
                )
            else:
                overpass = Overpass(
                    start_pt=start_pt,
                    max_pt=max_pt,
                    end_pt=end_pt
                )
            overpasses[j] = overpass


    def peakmem_predict_realtime_compute(self):
        location = self.location
        t = self.t
        sat = self.sat
        sun = self.sun

        store_sat_id = False
        
        min_elevation = 10.1

        rho = RhoVector()
        rho.time = t
        rsiteECEF = site_ECEF(location.lat, location.lon, location.h)
        rho.rECEF = site_sat_rotations(rsiteECEF, sat.rECEF)
        rho.rSEZ = ecef2sez(rho.rECEF, location.lat, location.lon)
        rng, az, el = razel(rho.rSEZ)
        rho.rng = rng
        rho.az = az
        rho.el = el
        # rsiteECI = ecef2eci(rsiteECEF, jdt_utc1)
        
        # Find Overpasses
        el0 = rho.el[:-1] - min_elevation
        el1 = rho.el[1:] - min_elevation
        el_change_sign = (el0*el1 < 0)   
        start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]  # Find the start of an overpass
        end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]    # Find the end of an overpass
        num_overpasses = min(start_idx.size, end_idx.size)       # Iterate over start/end indecies and gather inbetween indecies
        if start_idx.size < end_idx.size:
            end_idx = end_idx[1:]
        overpasses = [None] * num_overpasses
        for j in range(num_overpasses):
            # Store indecies of overpasses in a list
            idx0 = start_idx[j]
            idxf = end_idx[j]
            overpass_idx = np.arange(idx0, idxf+1, dtype=int)
            idxmax = np.argmax(el[overpass_idx])
            start_pt = Point(
                datetime=jday2datetime(rho.time.jd[idx0]),
                azimuth=rho.az[idx0],
                elevation=rho.el[idx0],
                range=rho.rng[idx0]
            )
            max_pt = Point(
                datetime=jday2datetime(rho.time.jd[idx0 + idxmax]),
                azimuth=rho.az[idx0 + idxmax],
                elevation=rho.el[idx0 + idxmax],
                range=rho.rng[idx0 + idxmax]
            )
            end_pt = Point(
                datetime=jday2datetime(rho.time.jd[idxf]),
                azimuth=rho.az[idxf],
                elevation=rho.el[idxf],
                range=rho.rng[idxf]
            )
            if store_sat_id:
                overpass = Overpass(
                    satellite_id=satellite.id,
                    start_pt=start_pt,
                    max_pt=max_pt,
                    end_pt=end_pt
                )
            else:
                overpass = Overpass(
                    start_pt=start_pt,
                    max_pt=max_pt,
                    end_pt=end_pt
                )
            overpasses[j] = overpass