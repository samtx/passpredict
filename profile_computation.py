# Profile Real-time computation of satellite overpasses

import cProfile
import pstats
import datetime
from pprint import pprint
import timeit

import numpy as np
import line_profiler

from passpredict.predictions import find_overpasses, compute_sun_data, compute_satellite_data
from passpredict.schemas import Location, Satellite, Tle
from passpredict.models import RhoVector, SatPredictData
from passpredict.timefn import compute_time_array_from_date
from passpredict.utils import epoch_from_tle


def f8(x):
    """ from https://gist.github.com/romuald/0346c76cfbbbceb3e4d1 """
    ret = "%8.3f" % x
    if ret != '   0.000':
        return ret
    return "%6dµs" % (x * int(1e6))
    
def f8_milli(x):
    ms = x * 1e3
    ret = f'{ms:6.2f}ms'
    return ret

def f8_micro(x):
    us = x * 1e6
    ret = f'{us:6.0f}µs'
    return ret


if __name__ == "__main__":


    # Set up satellite position
    dt_seconds = 1
    min_elevation = 10.0

    satellite = Satellite(id=25544, name='ISS')
    location = Location(lat=30.2711, lon=-97.7434, h=0, name='Austin, Texas')
    dt_seconds = 1.0
    tle1 = '1 25544U 98067A   20196.51422950 -.00000046  00000-0  72206-5 0  9999'
    tle2 = '2 25544  51.6443 213.2207 0001423 114.8006 342.8278 15.49514729236251'    
    tle = Tle(
        tle1=tle1,
        tle2=tle2,
        epoch=epoch_from_tle(tle1),
        satid=satellite.id
    ) 
    date_start = datetime.date(2020, 7, 14)
    date_end = date_start + datetime.timedelta(days=14)
    t = compute_time_array_from_date(date_start, date_end, dt_seconds)
    sun = compute_sun_data(t)
    sat = compute_satellite_data(tle, t, sun)
    sat = SatPredictData(id=sat.id, rECEF=sat.rECEF, illuminated=sat.illuminated)

    n = 50

    # benchmark
    timer = timeit.Timer(
        'find_overpasses(t, location, [sat], sun)',
        'from passpredict.predictions import find_overpasses',
        globals={
            'location': location,
            'sat': sat,
            't': t,
            'sun': sun
        }
    )
    try:
        res = timer.timeit(number=n)
    except Exception:
        timer.print_exc()
        res = 0

    avg_res = res/n
    if res > 0:
        print(f'\nAverage time: {avg_res*1e3:.2f}ms from {n} runs\n')
    
    if avg_res < 0.5:
        pstats.f8 = f8_milli # if benchmark is less than 1 second, return in microseconds
    else:
        pstats.f8 = f8
    
    # profile
    with cProfile.Profile() as pr:
        overpasses = find_overpasses(t, location, [sat], sun)

    pr.dump_stats('profile_computation.prof')
    
    ps = pstats.Stats(pr).sort_stats('cumulative')
    ps.print_stats()

    # Line Profile

    lp = line_profiler.LineProfiler()
    lp.add_function(RhoVector.find_overpasses)
    lp.runcall(find_overpasses, t, location, [sat], sun)
    lp_stats = lp.get_stats()
    lp.print_stats()
