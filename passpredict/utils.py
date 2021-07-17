from itertools import zip_longest
import shelve
import time

import numpy as np


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


class Cache:
    """
    A cache for downloaded TLE data using a shelf.
    """
    filename = 'tle.db'
        
    def __init__(self, filename=filename, ttl=86400):   
        self.filename = filename
        self.ttl = ttl
        self.cache = {}

    def set(self, key, value):
        key_hash = self._hash(key)
        self.cache[key_hash] = value
            
    def get(self, key, *a):
        key_hash = self._hash(key)
        value = self.cache.get(key_hash, *a)
        return value

    def pop(self, key, default_value=None, index=None):
        key_hash = self._hash(key)
        if key_hash in self.cache:
            value = self.get(key)
            del self.cache[key_hash]
            return value
        else:
            return default_value        

    def __contains__(self, key):
        key_hash = self._hash(key)
        return key_hash in self.cache

    def __setitem__(self, key, value):
        self.set(key, value)

    def __getitem__(self, key):
        return self.get(key)

    def __delitem__(self, key):
        key_hash = self._hash(key)
        del self.cache[key_hash]

    def _hash(self, key):
        return str(key)

    def __enter__(self):
        return self.open()

    def __exit__(self, *a):
        self.close()

    def open(self):
        self.cache = shelve.DbfilenameShelf(self.filename)
        return self
        
    def close(self):
        self.cache.close()
