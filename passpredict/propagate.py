
import numpy as np
import functools
from passpredict.timefn import julian_date, jdt_tsince, invjday, jday2datetime, \
                                 jday2npdatetime64, truncate_datetime
from passpredict.solar import is_sat_illuminated, sun_pos
from passpredict.models import SatelliteRV, Tle
from passpredict.rotations import teme2ecef
from passpredict.utils import parse_tles_from_celestrak, epoch_from_tle
import datetime
from sgp4.api import Satrec, WGS84
import json
import requests
import os

# use_cython = False
# try:
#     from passpredict._sgp4 import sgp4 as sgp4_pyx
#     use_cython = True
# except ImportError:
#     pass
# from passpredict.sgp4 import sgp4



def get_TLE(satellite, tle_data=None):
    if tle_data is None:
        if not os.path.exists('tle_data.json'):
            tle_data = parse_tles_from_celestrak()
            with open('tle_data.json', 'w') as file:
                json.dump(tle_data, file)
        else:
            with open('tle_data.json', 'r') as file:
                tle_data = json.load(file)
    tle1 = tle_data[str(satellite.id)]['tle1']
    tle2 = tle_data[str(satellite.id)]['tle2']
    epoch = epoch_from_tle(tle1)
    tle = Tle(tle1=tle1, tle2=tle2, epoch=epoch, satellite=satellite)
    return tle


def save_TLE_data(url=None):
    tle_data = parse_tles_from_celestrak(url)
    with open('tle_data.json', 'w') as file:
        json.dump(tle_data, file)


@functools.lru_cache(maxsize=64)
def propagate(tle1, tle2, dt0, dtf, dtsec=1.0):
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

    # # create array of julian dates to pass into sgp4
    jdt0 = julian_date(dt0)
    jdtf = julian_date(dtf)
    total_days = (dtf-dt0).total_seconds()/60
    dt_days = dtsec/(24*60*60.0)
    jdt = np.arange(jdt0, jdtf, dt_days, dtype=float)
    jd_array, fr_array = np.divmod(jdt, 1)

    error, r, v = satrec.sgp4_array(jd_array, fr_array)

    # Change arrays from column major to row major while keeping C-continuous
    r = np.reshape(r.ravel(order='F'), (3, r.shape[0]))
    v = np.reshape(v.ravel(order='F'), (3, v.shape[0]))

    # perform rotations to fixed earth coordinates
    rECI = r.copy()

    # Get earth orientation parameters deltaUTC1, xp, yp
    def eop(jdt):
        # Code this!!
        return 0, 0, 0

    deltaUTC1, xp, yp = eop(jdt)
    jdt_utc1 = jdt + deltaUTC1
    
    rECEF = teme2ecef(r, jdt_utc1, xp, yp)

    # Compute sun-satellite quantities
    # rsunECI = sun_pos(jdt)
    # sat_illum = is_sat_illuminated(rECI, rsunECI)
    dt_ary = jday2datetime(jdt)

    # Return output
    satellite_rv = SatelliteRV()
    satellite_rv.rECEF = rECEF
    satellite_rv.rECI = rECI
    # satellite_rv.visible = sat_illum
    satellite_rv.datetime = dt_ary
    satellite_rv.julian_date = jdt
    return satellite_rv


if __name__ == "__main__":
    # Test propagation routines
    from passpredict.models import Satellite

    # ISS satellite id
    satellite = Satellite(25544, "Int. Space Station")

    # ISS
    # tle1 = "1 25544U 98067A   19293.90487327  .00016717  00000-0  10270-3 0  9034"
    # tle2 = "2 25544  51.6426  97.8977 0006846 170.6875 189.4404 15.50212100 34757"

    tle = get_TLE(satellite)
    dt0 = datetime.datetime.now()#(2020, 5, 23, 0, 0, 0)
    dtf = dt0 + datetime.timedelta(days=14)

    print('start 1')
    satellite = propagate(tle.tle1, tle.tle2, dt0, dtf, dtsec=1)
    print('end 1')
    print('start 2')
    satellite = propagate(tle.tle1, tle.tle2, dt0, dtf, dtsec=10)
    print('end 2')

    # save position data
    r = satellite.rECI
    # np.save('r.npy', )

