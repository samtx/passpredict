import time
import datetime

import numpy as np

from passpredict import Location, Satellite, OMM
from passpredict.predict import compute_elevation_angle
from passpredict import predict_py
from passpredict.timefn import jday2datetime, julian_date
from passpredict.propagate import pkepler

# class CentralTime(datetime.tzinfo):
#     def __init__(self):
#         super().__init__(offset=datetime.timedelta(hours=-6))


def benchmark_compute_elevation_angle():
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    omm = OMM.from_tle(tle1, tle2)
    satellite = Satellite(omm)
    location = Location(32.1, -97.5, 20)
    jd = 2458871.5

    # measure execution time
    n = int(1e3)
    t1 = time.perf_counter()
    for _ in range(n):
        el = compute_elevation_angle(jd, location, satellite)
    t2 = time.perf_counter()
    print(f"Time: {(t2-t1)/n*1e6:f} us")


def predict_algorithm():
    """
    Test out main prediction algorithm
    """
    tle1 = '1 25544U 98067A   21258.74288194 -.00052450  00000-0 -97172-3 0  9998'
    tle2 = '2 25544  51.6428 254.9623 0002916  17.0077 340.1623 15.48363749302621'
    omm = OMM.from_tle(tle1, tle2)
    satellite = Satellite(omm)
    location = Location(30.1957, -97.8649, 0)
    observer = predict_py.Observer(location, satellite)
    t1 = time.perf_counter()
    res = observer.predict_overpasses()
    t2 = time.perf_counter()
    print(f'Time: {(t2-t1)*1000:f} ms')
    offset = datetime.timedelta(hours=-5)
    datetime_local = [d.replace(tzinfo=None) + offset for d in res]
    from pprint import pprint
    pprint(datetime_local)


def propagate_satellite():
    """
    Test out main prediction algorithm
    """
    tle1 = '1 25544U 98067A   21258.74288194 -.00052450  00000-0 -97172-3 0  9998'
    tle2 = '2 25544  51.6428 254.9623 0002916  17.0077 340.1623 15.48363749302621'
    satellite = Satellite.from_tle(tle1, tle2)
    jd0 = satellite.epoch
    jd = np.linspace(jd0, jd0 + 2, 1440 * 2)
    r = np.zeros((jd.size, 3))
    t1 = time.perf_counter()
    for i in range(jd.size):
        r[i] = satellite.propagate(jd[i])
    t2 = time.perf_counter()
    print(f'Time: {(t2-t1)*1000:f} ms')


def predict_example_11_6():
    """
    Using orbital elements for MIR space station, p. 106 Vallado
    """
    tle1 = "1 16609U 86017A   93352.53502934  .00007889  00000-0  10529-3 0   34"
    tle2 = "2 16609  51.6190  13.3340 0005770 102.5680 257.5950 15.5911407044786"
    satellite = Satellite.from_tle(tle1, tle2)
    r0 = np.array([6585.038266, 1568.184321, 9.116355])
    v0 = np.array([-1.1157766, 4.6316816, 6.0149576])
    epoch = 2450540.4
    satellite.set_epoch(epoch)
    jd = 2450540.5472
    dt_sec = (jd - epoch) * 86400
    r, v = pkepler(r0, v0, dt_sec, ndot=satellite.ndot, nddot=satellite.nddot)
    print(r)
    print(v)
    r, v = pkepler(r0, v0, dt_sec)
    print(r)
    print(v)
    print(satellite)

if __name__ == "__main__":
    # predict_algorithm()
    # propagate_satellite()
    # benchmark_compute_elevation_angle()
    predict_example_11_6()
