from passpredictor.sgp4io import twoline2rv, Satellite, wgs72, wgs84
from passpredictor.timefn import julian_date, jdt_tsince, invjday, jday2datetime, jday2npdatetime64
from passpredictor.rotations import TEME_to_ECEF
from passpredictor.predict import sun_pos, is_sat_illuminated, site_ECEF, razel
import passpredictor.rotations as rotations
import passpredictor.predict as predict
import numpy as np
import datetime

import matplotlib.pyplot as plt

use_cython = False
try:
    from passpredictor._sgp4 import sgp4 as sgp4_pyx
    use_cython = True
except ImportError:
    pass

from passpredictor.sgp4 import sgp4


ISS_TLE = (
"1 25544U 98067A   19293.90487327  .00016717  00000-0  10270-3 0  9034",
"2 25544  51.6426  97.8977 0006846 170.6875 189.4404 15.50212100 34757"
)

class Point(object):
    def __init__(self, datetime, azimuth, elevation, range_):
        self.datetime = datetime
        self.azimuth = azimuth
        self.elevation = elevation
        self.range = range_
    def __repr__(self):
        dtstr = self.datetime.strftime("%b %d %Y, %H:%M:%S")
        s = "{}UTC el={:.1f}d, az={:.1f}d, rng={:.1f}km".format(
            dtstr, self.elevation, self.azimuth, self.range)
        return s

class Overpass(object):
    def __init__(self, start_pt, max_pt, end_pt, t, r):
        self.start_pt = start_pt
        self.max_pt = max_pt
        self.end_pt = end_pt
        self.t = t
        self.r = r
        self.sat = None
        self.location = None

class Location(object):
    def __init__(self, lat, lon, h, name=None, tz=None):
        self.lat = lat
        self.lon = lon
        self.h = h
        self.name = name
        self.tz = tz
    def __repr__(self):
        return self.name

class SatelliteRV(object):
    def __init__(self):
        self.satellite = None
        self.tle = None
        self.rsun = None
        self.dt = None
        self.jdt = None
        self.rECEF = None
        self.rECI = None
        self.modified = None
        self.deltaT = None
        self.lat = None
        self.lon = None
        self.alt = None
        self.is_illum = None

class SunPosition(object):
    def __init__(self):
        self.rECI = None
        self.jdt = None
        self.dt_ary = None

class TLE(object):
    def __init__(self, tle1, tle2, dt):
        self.tle1 = tle1
        self.tle2 = tle2
        self.dt = dt
        self.satellite = None

class SatelliteDetail(object):
    def __init__(self, satid, name):
        self.id = satid
        self.name = name

def propagate_tle(satellite, sun):
    tle = satellite.tle
    satrec = twoline2rv(tle.tle1, tle.tle2, wgs84)
    min_per_day = 1440
    total_days = 10
    dt_per_min = 60
    t = np.linspace(0, min_per_day*total_days, min_per_day*total_days*dt_per_min)
    r = np.empty((3, t.size))
    v = np.empty((3, t.size))
    for i in range(t.size):
        if i % (min_per_day*dt_per_min) == 0:
            print(f'i = {i}')
        ri, vi = sgp4(satrec, t[i], wgs84)
        r[:, i] = ri
        v[:, i] = vi
    # print(satrec)
    epoch = julian_date(satrec.epoch)
    jdt = jdt_tsince(epoch, t)
    rECI = r
    rECEF = TEME_to_ECEF(r, jdt)

    # get array of corresponding datetime and np.datetime64 objects
    dt64_0 = jday2npdatetime64(jdt[0])
    dt64_step = jday2npdatetime64(jdt[1]) - dt64_0
    dt64_n = jday2npdatetime64(jdt[-1])
    dt64_ary = np.arange(dt64_0, dt64_n+dt64_step, dt64_step, dtype=np.datetime64)
    dt_ary = dt64_ary.astype(datetime.datetime)
    dt_step = dt64_step.astype(datetime.timedelta)

    # Compute sun-satellite quantities
    sat_illum = is_sat_illuminated(rECI, sun.rECI)

    np.save('dt_ary.npy', dt_ary, allow_pickle=True)

    np.savez('out.npz', rECEF=rECEF, rTEME=r, vTEME=v, tsince=t,
             tle=tle, jdt=jdt, epoch=epoch,
             dt64_ary=dt64_ary, dt64_step=dt64_step,
             dt_ary=dt_ary, dt_step=dt_step,
             rsun=rsun, sat_illum=sat_illum)


def site_rotations(lat, lon, h, rECEF):
    rsiteECEF = site_ECEF(lat, lon, h)
    rho = rECEF - np.atleast_2d(rsiteECEF).T
    rSEZ = rotations.ECEF_to_SEZ(rho, lat, lon)
    rng, az, el = razel(rSEZ)
    return rng, az, el


def realtime_compute(lat, lon, h, rECEF, jdt):
    rsiteECEF = predict.site_ECEF(lat, lon, h)
    rho = rECEF - np.atleast_2d(rsiteECEF).T
    rSEZ = rotations.ECEF_to_SEZ(rho, lat, lon)
    rng, az, el = predict.razel(rSEZ)
    np.savez('out2.npz', az=az, el=el, rng=rng, rSEZ=rSEZ, jdt=jdt)


def plot_razel():
    with np.load('out2.npz') as data:
        el = data['el']
        az = data['az']
        r = data['rSEZ']
        jdt = data['jdt']
    dtobjs = invjday(jdt)
    # dt =  invjday_to_npdatetime64(*dtobjs)
    # time zone shift UTC-7 for Tucson
    dt = np.datetime64(dtobjs)
    dt += datetime.timedelta(hours=-7)
    idx = r[2] > 0
    fig, ax = plt.subplots(2,1)
    ax[0].plot(dt, el, label='elevation')
    ax[0].legend()
    ax[0].grid()
    ax[1].plot(dt, az, label='azimuth')
    ax[1].legend()
    ax[1].grid()
    plt.show()

def get_overpasses(el, azm, rng, dt_ary, rhoSEZ, sat, loc):
    el0 = el[:-1]
    el1 = el[1:]
    el_change_sign = (el0*el1 < 0)
    # Find the start of an overpass
    start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]
    # Find the end of an overpass
    end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]
    # print(f'start shape {start_idx.shape}  end shape = start shape {end_idx.shape}')

    # Iterate over start/end indecies and gather inbetween indecies
    overpasses = np.empty(start_idx.size, dtype=object)
    for j in range(start_idx.size):
        # Store indecies of overpasses in a list
        idx0 = start_idx[j]
        idxf = end_idx[j]
        overpass_idx = np.arange(idx0, idxf+1, dtype=int)
        idxmax = np.argmax(el[overpass_idx])
        start_pt = Point(
            dt_ary[idx0],
            az[idx0],
            el[idx0],
            rng[idx0]
        ),
        max_pt = Point(
            dt_ary[idxmax],
            az[idxmax],
            el[idxmax],
            rng[idxmax]
        )
        end_pt = Point(
            dt_ary[idxf],
            az[idxf],
            el[idxf],
            rng[idxf]
        ),
        overpass = Overpass(
            start_pt,
            max_pt,
            end_pt,
            dt_ary[overpass_idx],
            rhoSEZ[:,overpass_idx]
        )
        overpass.location = loc
        overpasses[j] = overpass

    return overpasses


if __name__ == "__main__":

    # Tucson
    lat, lon, h = 32.2226, -110.9747, 0

    sat = Satellite()
    tle = TLE(ISS_TLE[0], ISS_TLE[1], None)
    sat.tle = tle
    print(tle.tle1, tle.tle2)

    jdt0 = np.array(2458777.40487327)
    sun = SunPosition()
    rsun = sun_pos(jdt0)


    # propagate_tle(sat)

    with np.load('out.npz') as data:
        rsatECEF = data['rECEF']
        jdt = data['jdt']

    dt_ary = np.load('dt_ary.npy', allow_pickle=True)

    loc = Location(lat, lon, h, 'Tucson')
    # realtime_compute(loc.lat, loc.lon, loc.h, rsatECEF, jdt)

    with np.load('out2.npz') as data:
        el = data['el']
        az = data['az']
        rng = data['rng']
        rSEZ = data['rSEZ']

    # overpasses = get_overpasses(el, az, rng, dt_ary, rSEZ, sat, loc)
    # np.save('overpasses.npy', overpasses, allow_pickle=True)

    overpasses = np.load('overpasses.npy', allow_pickle=True)
    print('hello')
    # realtime_compute(lat, lon, h, rECEF, jdt)
    # plot_razel()