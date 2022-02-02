import json
import shelve
import time
import pathlib
import abc


class BaseCache(abc.ABC):
    """
    Base object for Caches
    """

    def _hash(self, key):
        return str(key)

    def __contains__(self, key):
        key_hash = self._hash(key)
        return key_hash in self.cache

    def get(self, key, ignore_ttl=False):
        if key not in self: return None
        item = self._get(key)
        if (not ignore_ttl) and (ttl := item['ttl']):
            if time.time() > ttl:
                self._del(key)
                return None
        return item['data']

    def set(self, key, value, ttl=None):
        if ttl is not None:
            ttl += time.time()
        item = {'ttl': ttl, 'data': value}
        self._set(key, item)

    @abc.abstractmethod
    def _set(self, key, value):
        raise NotImplementedError

    @abc.abstractmethod
    def _get(self, key):
        raise NotImplementedError

    @abc.abstractmethod
    def _del(self, key):
        raise NotImplementedError

    def pop(self, key, default_value=None):
        key_hash = self._hash(key)
        if key_hash in self.cache:
            value = self.get(key)
            self._del(key_hash)
            return value
        else:
            return default_value

    def __enter__(self):
        self.load()
        return self

    def __exit__(self, *a):
        self.save()

    @abc.abstractmethod
    def load(self):
        raise NotImplementedError

    @abc.abstractmethod
    def save(self):
        raise NotImplementedError


class JsonCache(BaseCache):
    """
    A cache for downloaded TLE data using a json file
    """
    filename = 'passpredict.json'

    def __init__(self, filename=filename):
        self.filename = filename
        self.cache = {}

    def _set(self, key, value):
        self.cache[key] = value

    def _get(self, key):
        return self.cache[key]

    def _del(self, key):
        del self.cache[key]

    def load(self, strict=False):
        """ Load cache from json file """
        try:
            with open(self.filename, 'r') as f:
                self.cache = json.load(f)
        except FileNotFoundError as e:
            if strict:
                raise e
            else:
                pass

    def save(self):
        """ Save cache to json file """
        with open(self.filename, 'w') as f:
            json.dump(self.cache, f, indent=2)


class MemoryCache(BaseCache):
    """
    A dictionary in-memory cache
    """
    def __init__(self):
        self.cache = {}

    def _set(self, key, value):
        self.cache[key] = value

    def _get(self, key):
        return self.cache[key]

    def _del(self, key):
        del self.cache[key]

    def load(self):
        pass

    def save(self):
        pass


class ShelfCache:
    """
    A cache for downloaded TLE data using a shelf.
    """
    filename = 'passpredict.db'

    def __init__(self, filename=filename):
        self.filename = filename

    def set(self, key, value):
        key_hash = self._hash(key)
        self.cache[key_hash] = value

    def get(self, key, *a):
        key_hash = self._hash(key)
        value = self.cache.get(key_hash, *a)
        return value

    def pop(self, key, default_value=None):
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
