# Profile Real-time computation of satellite overpasses

import cProfile
import datetime
from pprint import pprint

import numpy as np

from passpredict.predictions import predict_passes
from passpredict.propagate import propagate
from passpredict.models import Location, Satellite
from passpredict.timefn import truncate_datetime
from passpredict.utils import get_TLE

def main():

    # Set up satellite position
    dt_seconds = 1
    num_days = 21
    min_elevation = 10.0

    austin = Location(lat=30.2672, lon=-97.7431, h=0, name='Austin')
    iss = Satellite(id=25544, name='ISS')
    iss_tle = get_TLE(iss)
    datetime_start = truncate_datetime(
        datetime.datetime.now(tz=datetime.timezone.utc)
    )
    datetime_end = datetime_start + datetime.timedelta(days=num_days)

    iss_rv = propagate.__wrapped__(
        iss_tle.tle1, iss_tle.tle2, datetime_start, datetime_end, dt_seconds
    )

    with cProfile.Profile() as pr:
        overpasses = predict_passes(
            austin.lat, austin.lon, austin.h,
            iss_rv.rECEF, iss_rv.rECI, iss_rv.julian_date,
            min_elevation=min_elevation#, loc=location, sat=satellite
        )

    pr.print_stats()
    pr.dump_stats('profile_computation.prof')


def main_lineprofile():
    # Set up satellite position
    dt_seconds = 1
    num_days = 14
    min_elevation = 10.0

    austin = Location(lat=30.2672, lon=-97.7431, h=0, name='Austin')
    iss = Satellite(id=25544, name='ISS')
    iss_tle = get_TLE(iss)
    datetime_start = truncate_datetime(
        datetime.datetime.now(tz=datetime.timezone.utc)
    )
    datetime_end = datetime_start + datetime.timedelta(days=num_days)

    iss_rv = propagate.__wrapped__(
        iss_tle.tle1, iss_tle.tle2, datetime_start, datetime_end, dt_seconds
    )

    overpasses = predict_passes(
        austin.lat, austin.lon, austin.h,
        iss_rv.rECEF, iss_rv.rECI, iss_rv.julian_date,
        min_elevation=min_elevation#, loc=location, sat=satellite
    )
    return 1


if __name__=="__main__":
    main_lineprofile()
    # main()
