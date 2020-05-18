
import numpy as np
import functools
from passpredictor.timefn import julian_date, jdt_tsince, invjday, jday2datetime, \
                                 jday2npdatetime64, truncate_datetime
from passpredictor.solar import is_sat_illuminated, sun_pos
from passpredictor.models import SatelliteRV, Tle
from passpredictor.rotations import teme2ecef
import datetime
from sgp4.api import Satrec, WGS84
from itertools import zip_longest
import json
import requests
import os

# use_cython = False
# try:
#     from passpredictor._sgp4 import sgp4 as sgp4_pyx
#     use_cython = True
# except ImportError:
#     pass
# from passpredictor.sgp4 import sgp4


def grouper(iterable, n, fillvalue=None):
    """
    from itertools recipes https://docs.python.org/3.7/library/itertools.html#itertools-recipes
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def get_TLE(satellite):
    if not os.path.exists('tle_data.json'):
        tle_data = save_TLE_data()
    else:
        with open('tle_data.json', 'r') as file:
            tle_data = json.load(file)
    tle1 = tle_data[str(satellite.id)]['tle1']
    tle2 = tle_data[str(satellite.id)]['tle2']
    epoch_year = int(tle1[18:20])
    if epoch_year < 57:
        epoch_year += 2000
    else:
        epoch_year += 1900
    epoch_day = float(tle1[20:32])
    epoch_day, epoch_day_fraction = np.divmod(epoch_day, 1)
    epoch_microseconds = epoch_day_fraction * 24 * 60 * 60 * 1e6
    epoch = datetime.datetime(epoch_year, month=1, day=1) + \
            datetime.timedelta(days=int(epoch_day-1)) + \
            datetime.timedelta(microseconds=int(epoch_microseconds))
    tle = Tle(tle1, tle2, epoch, satellite)
    return tle


def save_TLE_data():
    """
    Download current TLEs from Celestrak and save them to a JSON file
    """
    stations_url = 'https://celestrak.com/NORAD/elements/stations.txt'
    tle_data = {}
    r = requests.get(stations_url, stream=True)
    for lines in grouper(r.iter_lines(), 3):
        tle0, tle1, tle2 = [line.decode('ascii') for line in lines]
        name = tle0.strip()  # satellite name
        satellite_id = tle1[2:7]
        tle_data.update(
            {
                satellite_id : {'name': name, 'tle1': tle1, 'tle2': tle2}
            }
        )
    with open('tle_data.json', 'w') as file:
        json.dump(tle_data, file)
    return tle_data


@functools.lru_cache(maxsize=64)
def propagate(tle1, tle2, dt0=None, dtf=None, dtsec=1.0):
    """Propagate satellite position forward in time.

    Parameters:
        tle1 : str
            first line of two line element set
        tle2 : str
            second line of two line element set
        dt0 : datetime
            initial datetime to begin propagation, to nearest second
        dtf : datetime
            final datetime for propagation, to nearest second
        dtsec : float
            time interval in seconds for propagation

    Returns:
        r : float (3, n)
            satellite position vector in TEME coordinates
        v : float (3, n)
            satellite velocity vector in TEME coordinates
    """

    # if not use_cython:
    #     sgp4fn = sgp4
    # else:
    #     sgp4fn = sgp4_pyx

    # Get TLE data for ISS
    # tle1, tle2 = get_TLE()

    satrec = Satrec.twoline2rv(tle1, tle2, WGS84)

    # truncate start and end time to the second
    # if dt0 is None:
    #     dt0 = satrec.epoch
    # if dtf is None:
    #     dtf = dt0 + datetime.timedelta(days=5)
    # dt0 = truncate_datetime(dt0)
    # dtf = truncate_datetime(dtf)

    # # find total minutes between start and end times
    # totmin = (dtf - dt0).total_seconds()/60.0

    # # create array of julian dates to pass into sgp4
    # jdt0 = satrec.jdsatepoch
    jdt0 = julian_date(dt0)
    jdtf = julian_date(dtf)
    total_days = (dtf-dt0).total_seconds()/60

    # t0 = (jdt0 - satrec.jdtsatepoch).total_seconds()/60.0
    dt_days = dtsec/(24*60*60.0)
    jdt = np.arange(jdt0, jdtf, dt_days, dtype=float)
    # dt_array =
    jd_array, fr_array = np.divmod(jdt, 1)

    # r = np.empty((3, t.size))
    # v = np.empty((3, t.size))
    # for i in range(t.size):
    #     # if i % (min_per_day*dt_per_min) == 0:
    #     #     print(f'i = {i}')
    #     ri, vi = sgp4fn(satrec, t[i], wgs84)
    #     r[:, i] = ri
    #     v[:, i] = vi

    # breakpoint()
    # r = np.empty((3, jd_array.size) , order='F')
    # v = np.empty((3, jd_array.size) , order='F')
    error, r, v = satrec.sgp4_array(jd_array, fr_array)

    # Change arrays from column major to row major while keeping C-continuous
    r = np.reshape(r.ravel(order='F'), (3, r.shape[0]))
    v = np.reshape(v.ravel(order='F'), (3, v.shape[0]))

    # perform rotations to fixed earth coordinates
    # jdt = jdt_tsince(jdt0, t)
    rECI = r.copy()
    rECEF = teme2ecef(r, jdt)

    # get array of corresponding datetime and np.datetime64 objects
    # dt64_0 = jday2npdatetime64(jdt[0])
    # dt64_step = jday2npdatetime64(jdt[1]) - dt64_0
    # dt64_n = jday2npdatetime64(jdt[-1])
    # dt64_ary = np.arange(dt64_0, dt64_n+dt64_step, dt64_step, dtype=np.datetime64)
    # dt_ary = dt64_ary.astype(datetime.datetime)
    # dt_step = dt64_step.astype(datetime.timedelta)

    # Compute sun-satellite quantities
    rsunECI = sun_pos(jdt)
    sat_illum = is_sat_illuminated(rECI, rsunECI)
    dt_ary = jday2datetime(jdt)

    # Return output
    satellite_rv = SatelliteRV()
    satellite_rv.rECEF = rECEF
    satellite_rv.rECI = rECI
    satellite_rv.visible = sat_illum
    satellite_rv.datetime = dt_ary
    # satellite.dt64 = dt64_ary
    satellite_rv.julian_date = jdt

    # sun = Sun()
    # sun.rECI = rsunECI
    # sun.jdt = jdt
    # sun.dt = dt_ary
    # sun.dt64 = dt64_ary

    return satellite_rv

    # np.save('dt_ary.npy', dt_ary, allow_pickle=True)
    # np.savez('out.npz', rECEF=rECEF, rTEME=r, vTEME=v, tsince=t,
    #          tle=tle, jdt=jdt, epoch=epoch,
    #          dt64_ary=dt64_ary, dt64_step=dt64_step,
    #          dt_ary=dt_ary, dt_step=dt_step,
    #          rsun=rsun, sat_illum=sat_illum)

if __name__ == "__main__":
    # Test propagation routines
    from passpredictor.models import Satellite
    save_TLE_data()

    # ISS satellite id
    satellite = Satellite(25544, "Int. Space Station")

    # ISS
    # tle1 = "1 25544U 98067A   19293.90487327  .00016717  00000-0  10270-3 0  9034"
    # tle2 = "2 25544  51.6426  97.8977 0006846 170.6875 189.4404 15.50212100 34757"

    tle = get_TLE(satellite)

    print('start 1')
    satellite = propagate(tle.tle1, tle.tle2, dtsec=1)
    print('end 1')
    print('start 2')
    satellite = propagate(tle.tle1, tle.tle2, dtsec=10)
    print('end 2')

    # save position data
    r = satellite.rECI
    # np.save('r.npy', )

