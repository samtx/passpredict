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