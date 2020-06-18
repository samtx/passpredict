import json
import datetime
from itertools import zip_longest
from collections import OrderedDict

import numpy as np

CACHE_DIRECTORY = ".passpredict_cache"


def grouper(iterable, n, fillvalue=None):
    """
    from itertools recipes https://docs.python.org/3.7/library/itertools.html#itertools-recipes
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def epoch_from_tle(tle1):
    """
    Extract epoch as datetime from tle line 1
    """
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
    return epoch


def get_orbit_data_from_celestrak(satellite_id):
    """

    Params:
        satellite_id : int
            NORAD satellite ID


    See https://celestrak.com/NORAD/documentation/gp-data-formats.php

    Can use the new celestrak api for satellite ID
    https://celestrak.com/NORAD/elements/gp.php?CATNR=25544&FORMAT=json

    for tle api:
    https://celestrak.com/satcat/tle.php?CATNR=25544

    Supplemental TLEs available: (not fully working as json)
    https://celestrak.com/NORAD/elements/supplemental/gp-index.php?GROUP=iss&FORMAT=json

    https://celestrak.com/NORAD/elements/supplemental/starlink.txt
    https://celestrak.com/NORAD/elements/supplemental/iss.txt
    
    """
    query = {
        'CATNR': satellite_id,
        'FORMAT': 'json'
    }
    url = 'https://celestrak.com/NORAD/elements/gp.php'
    r = requests.get(url, data=query)
    return r.json()


def parse_tles_from_celestrak(url=None):
    """
    Download current TLEs from Celestrak and save them to a JSON file
    
    """
    if url is None:
        url = 'https://celestrak.com/NORAD/elements/stations.txt'
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
    return tle_data


class Cache():
    def __init__(self, cache_directory=CACHE_DIRECTORY, ttl=30):
        self.directory = cache_directory
        self.ttl = ttl # days
        self.cache = OrderedDict()

    def set(self, key, val):
        pass

    def get(self, key, val):
        return None


def cache_tle(data, method, cache_directory=CACHE_DIRECTORY):
    """
    Get/Set satellite TLE data

    Params:
        data = tle data
        method = str ['get','set']
        cache_directory: str, optional
    """
    pass


def cache_satellite_position(data, method, cache_directory=CACHE_DIRECTORY):
    """
    Get/Set satellite ECEF position vectors in cache
    """
    pass
