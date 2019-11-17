
import numpy as np
import functools
from passpredictor.sgp4io import twoline2rv, Satellite, wgs72, wgs84
from passpredictor.timefn import julian_date, jdt_tsince, invjday, jday2datetime, \
                                 jday2npdatetime64, truncate_datetime
from passpredictor.predict import is_sat_illuminated, sun_pos
from passpredictor.models import SatelliteRV, SunPosition
from passpredictor.rotations import teme2ecef
import datetime


use_cython = False
try:
    from passpredictor._sgp4 import sgp4 as sgp4_pyx
    use_cython = True
except ImportError:
    pass
from passpredictor.sgp4 import sgp4

@functools.lru_cache(maxsize=64)
def propagate(tle1, tle2, dt0=None, dtf=None, dtsec=1.0, use_cython=use_cython):
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

    if not use_cython:
        sgp4fn = sgp4
    else:
        sgp4fn = sgp4_pyx

    satrec = twoline2rv(tle1, tle2, wgs84)

    # get start and end time to neareast second

    # truncate start and end time to the second
    if dt0 is None:
        dt0 = satrec.epoch
    if dtf is None:
        dtf = dt0 + datetime.timedelta(days=5)
    dt0 = truncate_datetime(dt0)
    dtf = truncate_datetime(dtf)

    # find total minutes between start and end times
    totmin = (dtf - dt0).total_seconds()/60.0

    # create array t of propagation minutes
    t0 = (dt0 - satrec.epoch).total_seconds()/60.0
    dtmin = dtsec/60.0
    tf = t0 + totmin + dtmin
    t = np.arange(t0, tf, dtmin, dtype=float)

    r = np.empty((3, t.size))
    v = np.empty((3, t.size))
    for i in range(t.size):
        # if i % (min_per_day*dt_per_min) == 0:
        #     print(f'i = {i}')
        ri, vi = sgp4fn(satrec, t[i], wgs84)
        r[:, i] = ri
        v[:, i] = vi

    # perform rotations to fixed earth coordinates
    jdt0 = julian_date(dt0)
    jdt = jdt_tsince(jdt0, t)
    rECI = r.copy()
    rECEF = teme2ecef(r, jdt)

    # get array of corresponding datetime and np.datetime64 objects
    dt64_0 = jday2npdatetime64(jdt[0])
    dt64_step = jday2npdatetime64(jdt[1]) - dt64_0
    dt64_n = jday2npdatetime64(jdt[-1])
    dt64_ary = np.arange(dt64_0, dt64_n+dt64_step, dt64_step, dtype=np.datetime64)
    dt_ary = dt64_ary.astype(datetime.datetime)
    dt_step = dt64_step.astype(datetime.timedelta)

    # Compute sun-satellite quantities
    rsunECI = sun_pos(jdt)
    sat_illum = is_sat_illuminated(rECI, rsunECI)

    # Return output
    sat = SatelliteRV()
    sat.tle1 = tle1
    sat.tle2 = tle2
    sat.rECEF = rECEF
    sat.rECI = rECI
    sat.is_illum = sat_illum
    sat.dt = dt_ary
    sat.dt64 = dt64_ary
    sat.jdt = jdt

    sun = SunPosition()
    sun.rECI = rsunECI
    sun.jdt = jdt
    sun.dt = dt_ary
    sun.dt64 = dt64_ary

    out = {
        'sat': sat,
        'sun': sun,
    }
    return out

    # np.save('dt_ary.npy', dt_ary, allow_pickle=True)
    # np.savez('out.npz', rECEF=rECEF, rTEME=r, vTEME=v, tsince=t,
    #          tle=tle, jdt=jdt, epoch=epoch,
    #          dt64_ary=dt64_ary, dt64_step=dt64_step,
    #          dt_ary=dt_ary, dt_step=dt_step,
    #          rsun=rsun, sat_illum=sat_illum)

if __name__ == "__main__":
    # Test propagation routines

    # ISS
    tle1 = "1 25544U 98067A   19293.90487327  .00016717  00000-0  10270-3 0  9034"
    tle2 = "2 25544  51.6426  97.8977 0006846 170.6875 189.4404 15.50212100 34757"

    print('start 1')
    out = propagate(tle1, tle2, dtsec=10)
    sat = out['sat']
    sun = out['sun']
    assert np.all(sat.dt == sun.dt)
    assert np.all(sat.jdt == sun.jdt)
    print('end 1')
    print('start 2')
    out = propagate(tle1, tle2, dtsec=10)
    print('end 2')

