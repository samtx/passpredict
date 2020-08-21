import json
import datetime
from itertools import zip_longest
from typing import NamedTuple
import shelve
import time
from pathlib import Path
# from collections.abc import Mapping

import numpy as np
import requests

from .schemas import Tle


def shift_angle(x: float) -> float:
    """Shift angle in radians to [-pi, pi)
    
    Args:
        x: float, angle in radians

    Reference: 
        https://stackoverflow.com/questions/15927755/opposite-of-numpy-unwrap/32266181#32266181
    """
    return (x + np.pi) % (2 * np.pi) - np.pi
    


def grouper(iterable, n, fillvalue=None):
    """
    from itertools recipes https://docs.python.org/3.7/library/itertools.html#itertools-recipes
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def epoch_from_tle_datetime(epoch_string: str) -> datetime.datetime:
    """
    Return datetime object from tle epoch string
    """
    epoch_year = int(epoch_string[0:2])
    if epoch_year < 57:
        epoch_year += 2000
    else:
        epoch_year += 1900
    epoch_day = float(epoch_string[2:])
    epoch_day, epoch_day_fraction = np.divmod(epoch_day, 1)
    epoch_microseconds = epoch_day_fraction * 24 * 60 * 60 * 1e6
    epoch = datetime.datetime(epoch_year, month=1, day=1) + \
            datetime.timedelta(days=int(epoch_day-1)) + \
            datetime.timedelta(microseconds=int(epoch_microseconds))
    return epoch
    

def epoch_from_tle(tle1: str) -> datetime.datetime:
    """
    Extract epoch as datetime from tle line 1
    """
    epoch_string = tle1[18:32]
    return epoch_from_tle_datetime(epoch_string)
    

def satid_from_tle(tle1: str) -> int:
    """
    Extract satellite NORAD ID as int from tle line 1
    """
    return int(tle1[2:7])


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


def parse_tles_from_celestrak(satellite_id=None):
    """
    Download current TLEs from Celestrak and save them to a JSON file
    
    """
    if satellite_id is None:
        url = 'https://celestrak.com/NORAD/elements/stations.txt'
        params = {}
    else:
        url = 'https://celestrak.com/satcat/tle.php'
        params = {'CATNR': satellite_id}
    r = requests.get(url, params=params, stream=True)
    tle_data = {}
    for tle_strings in grouper(r.text.splitlines(), 3):
        tle_data.update(parse_tle(tle_strings))
    return tle_data


def parse_tle(tle_string_list):
    """
    Parse a single 3-line TLE from celestrak
    """
    tle0, tle1, tle2 = tle_string_list
    name = tle0.strip()  # satellite name
    satellite_id = satid_from_tle(tle1)
    return {satellite_id : {'name': name, 'tle1': tle1, 'tle2': tle2}}


def get_TLE(satid: int, tle_data=None):
    tle_data = parse_tles_from_celestrak(satid)
    tle1 = tle_data[satid]['tle1']
    tle2 = tle_data[satid]['tle2']
    epoch = epoch_from_tle(tle1)
    tle = Tle(tle1=tle1, tle2=tle2, epoch=epoch, satid=satid)
    return tle


def save_TLE_data(url=None):
    tle_data = parse_tles_from_celestrak(url)
    with open('tle_data.json', 'w') as file:
        json.dump(tle_data, file)


class CacheItem(NamedTuple):
    data: object
    ttl: int = -1    # time to live in seconds


class Cache:
    cache_filename = 'passpredict_cache.db'
        
    def __init__(self, filename='passpredict_cache.db', ttl=84600):   
        self.filename = filename
        self.ttl_default = ttl
        self.cache = {}
        self.category = None

    def _get_ttl_timestamp(self, ttl: int) -> int:
        return int(time.time() + ttl)

    def set(self, key, value, ttl: int = None):
        key_hash = self.hash(key)
        if ttl is None:
            ttl = self.ttl_default
        ttl_timestamp = self._get_ttl_timestamp(ttl)
        self.cache[key_hash] = CacheItem(data=value, ttl=ttl_timestamp)
            
    def get(self, key, *a):
        ttl_now = time.time()
        key_hash = self.hash(key)
        item = self.cache.get(key_hash, *a)
        if (item is None) or (0 < item.ttl < ttl_now):
            return None
        else:
            return item.data

    def pop(self, key, default_value=None):
        key_hash = self.hash(key)
        if key_hash in self.cache:
            value = self.get(key)
            del self.cache[key_hash]
            return value
        else:
            return default_value

    def __contains__(self, key):
        key_hash = self.hash(key)
        return key_hash in self.cache

    def __setitem__(self, key, value):
        self.set(key, value)

    def __getitem__(self, key):
        return self.get(key)

    def __delitem__(self, key):
        key_hash = self.hash(key)
        del self.cache[key_hash]

    def hash(self, key):
        return str(key)

    def flush(self):
        """Remove expired cache entires"""
        ttl_now = time.time()
        for key, item in list(self.cache.items()):
            if 0 < item.ttl < ttl_now:
                key_hash = self.hash(key)
                del self.cache[key_hash]


    def __enter__(self):
        return self.open()

    def __exit__(self, *a):
        self.close()

    def open(self):
        """
        Need to create the cache directory if it doesn't exist
        This is due to a bug that is fixed in PR 20274 but isn't merged yet [https://github.com/python/cpython/pull/20274]
        """
        filename_path = Path(self.filename)
        dir_path = filename_path.parent
        dir_path.mkdir(parents=True, exist_ok=True)
        self.cache = shelve.DbfilenameShelf(str(filename_path))
        return self
        
    def close(self):
        self.cache.close()