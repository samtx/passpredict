
import numpy as np
from passpredictor.sgp4 import cpropagation

cython_installed = False
try:
    from sgp4.cpropagation import sgp4
except ImportError:
    from sgp4.propagation import sgp4


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