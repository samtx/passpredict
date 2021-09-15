import time

from passpredict import Location, Satellite, OMM
from passpredict.predict import compute_elevation_angle

def main():
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


if __name__ == "__main__":
    main()
