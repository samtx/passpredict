from passpredict.sgp4io import twoline2rv, Satellite, wgs72, wgs84
from passpredict.timefn import julian_date, jdt_tsince, invjday, jday2datetime, jday2npdatetime64
from passpredict.rotations import teme2ecef
from passpredict.predict import sun_pos, is_sat_illuminated, site_ECEF, razel
import passpredict.rotations as rotations
import passpredict.predict as predict
from passpredict.propagate import propagate
from passpredict.models import SatelliteRV, SunPosition, Tle
import numpy as np
import datetime

import matplotlib.pyplot as plt

use_cython = False
try:
    from passpredict._sgp4 import sgp4 as sgp4_pyx
    use_cython = True
except ImportError:
    pass

from passpredict.sgp4 import sgp4


def realtime_compute(lat, lon, h, rECEF, jdt):
    rsiteECEF = predict.site_ECEF(lat, lon, h)
    rho = rECEF - np.atleast_2d(rsiteECEF).T
    rSEZ = rotations.ecef2sez(rho, lat, lon)
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


if __name__ == "__main__":

    # Tucson
    lat, lon, h = 32.2226, -110.9747, 0


    iss_tle = (
    "1 25544U 98067A   19293.90487327  .00016717  00000-0  10270-3 0  9034",
    "2 25544  51.6426  97.8977 0006846 170.6875 189.4404 15.50212100 34757"
    )

    tle = Tle(iss_tle[0], iss_tle[1], None)

    sat, sun = propagate(tle.tle1, tle.tle2, dtsec=10)
    sat.tle = tle

    loc = Location(lat, lon, h, 'Tucson')
    predict.realtime_compute(loc.lat, loc.lon, loc.h, rsatECEF, jdt)

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