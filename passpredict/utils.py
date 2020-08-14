import json
import datetime
import os
from itertools import zip_longest
from collections import OrderedDict
from typing import NamedTuple
import shelve
import time
from pathlib import Path
from enum import Enum
try:
    import cPickle as pickle
except:
    import pickle

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
    satellite_id = tle1[2:7]
    return {satellite_id : {'name': name, 'tle1': tle1, 'tle2': tle2}}


def get_TLE(satellite, tle_data=None):
    tle_data = parse_tles_from_celestrak(satellite.id)
    tle1 = tle_data[str(satellite.id)]['tle1']
    tle2 = tle_data[str(satellite.id)]['tle2']
    epoch = epoch_from_tle(tle1)
    tle = Tle(tle1=tle1, tle2=tle2, epoch=epoch, satellite=satellite)
    return tle


def save_TLE_data(url=None):
    tle_data = parse_tles_from_celestrak(url)
    with open('tle_data.json', 'w') as file:
        json.dump(tle_data, file)


class ShelfCache(shelve.DbfilenameShelf):
    """
    Default cache backend 
    """

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)

    def set(self, key, value):
        self.__setitem__(key, value)

    # def open(self, *a, **kw):
    #     self.cache = shelve.open(*a, **kw)


class DataCategory(str, Enum):
    time = 'time'
    sun = 'sun'
    tle = 'tle'
    sat = 'sat'


class CacheItem(NamedTuple):
    data: object
    ttl: int = None            # time to live in seconds, default is one day
    timestamp: int = None   # unix timestamp
    category: DataCategory = None


class Cache:
    data_cache_filename = 'cache_data.db'
    indices = {  # index table name: index item name
        'ttl': 'timestamp',
        'cat': 'category',
      }
        
    def __init__(self, cache_directory='.passpredict_cache', ttl=84600):   
        self.directory = cache_directory
        self.dir_path = Path(cache_directory)
        self.ttl_default = ttl
        self.cache = None
        self.category = None

    def _get_ttl_timestamp(self, ttl: int) -> int:
        return int(time.time() + ttl)

    def set(self, key, value, ttl: int = 0, category: DataCategory = None):
        key_hash = self.hash(key)
        timestamp = self._get_ttl_timestamp(ttl)
        # if key_hash in self.cache:
        #     old_item = self.cache[key_hash]
        #     # Remove old item from indices
        #     for index_tbl, index_field in self.indices:
        #         self.pop_index_item(index_tbl, old_item.__dict__() )
        #     old_timestamp = old_item.timestamp
        #     if old_timestamp is not None:
        #         self.ttl[old_timestamp].remove(key_hash)
        self.cache[key_hash] = CacheItem(data=value, ttl=ttl, timestamp=timestamp, category=category)
        # if ttl > 0:
        #     self._set_ttl(key_hash, timestamp)
        
            
    def _set_ttl(self, key_hash, timestamp: int):
        timestamp_keys = self.ttl.get(timestamp)
        if timestamp_keys is not None:
            timestamp_keys.add(key_hash)
            self.ttl[timestamp] = timestamp_keys
        else:
            self.ttl[timestamp] = {key_hash}

    def get(self, key, *a):
        key_hash = self.hash(key)
        item = self.cache.get(key_hash, *a)
        if (item is None) or ((item.ttl is not None) and (time.time() > item.timestamp)):
            return None
        else:
            return item.data

    def pop(self, key, default_value=None, index=None):
        key_hash = self.hash(key)
        if key_hash in self.cache:
            value = self.get(key)
            del self.cache[key_hash]
            # if index is not None:
            #     self.pop_index_item(index, )
            return value
        else:
            return default_value

    def get_index_items(self, index, index_key):
        index_dict = self.cache[index + '_index']
        return index_dict.get(index_key)

    def set_index_item(self, index, index_key, index_item):
        key = index + '_index'
        index_dict = self.cache[key]
        index_result = index_dict.get(index_key, set())
        index_result += {index_item}
        index_dict.update({
            index_key: index_result
        })
        self.cache[key] = index_dict

    def remove_item_from_cache(self, item):
        pass

    def pop_index_item(self, index, index_key, index_item):
        key = index + '_index'
        index_dict = self.cache[key]
        index_result = index_dict.get(index_key, set())
        if index_item in index_result:
            index_result.remove(index_item)
        if len(index_result) == 0:
            index_dict.pop(index_key)
        else:
            index_dict.update({
                index_key: index_result
            })
        self.cache[key] = index_dict

        

    def __contains__(self, key):
        key_hash = self.hash(key)
        return key_hash in self.cache

    def __setitem__(self, key, value):
        self.set(key, value)

    def __getitem__(self, key):
        return self.get(key)

    def hash(self, key):
        return str(key)

    def flush(self):
        """Remove expired cache entires"""
        current_time = time.time()
        for key, value in self.cache.items():
            t = value.timestamp
            if (t is not None) and (current_time > t):
                del self.cache[key]
        # ttl_index_dict = self.cache['ttl_index']
        # current_timestamp = time.time()
        # expired_keys = (key for key in ttl_index_dict.keys() if key >= current_timestamp)
        # for key in expired_keys:


    def __enter__(self, *a, **kw):
        return self.open()

    def __exit__(self, *a, **kw):
        self.close()

    def _build_indices(self):
        for idx, _ in self.indices:
            key = idx + '_index'
            self.cache[key] = OrderedDict()

    # def _open_ttl(self):
    #     self.ttl_path = self.dir_path / self.data_index_filename
    #     if self.ttl_path.exists():
    #         with open(self.ttl_path, 'rb+') as f:
    #             self.ttl = pickle.load(f)
    #     else:
    #         self.ttl = OrderedDict()
    #         self._build_ttl_index()

    def open(self):
        """
        Need to create the cache directory if it doesn't exist
        This is due to a bug that is fixed in PR 20274 but isn't merged yet [https://github.com/python/cpython/pull/20274]
        """
        dir_path = Path(self.directory)
        if not dir_path.exists():
            dir_path.mkdir()
        cache_path = dir_path / self.data_cache_filename
        self.cache = shelve.DbfilenameShelf(str(cache_path))
        # self._open_ttl()
        return self
        
    def close(self):
        self.cache.close()
        # with open(self.ttl_path, 'wb') as f:
        #     pickle.dump(self.ttl, f)

    def _build_ttl_index(self):
        """Iterate through cache and build ttl index"""
        assert self.ttl is not None
        pass
