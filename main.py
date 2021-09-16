import time
import datetime

from passpredict import Location, Satellite, OMM
from passpredict.predict import compute_elevation_angle
from passpredict import predict_py

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


if __name__ == "__main__":
    predict_algorithm()
    # benchmark_compute_elevation_angle()
