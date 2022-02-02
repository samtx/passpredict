import json
import shelve
import time
import abc


class BaseCache(abc.ABC):
    """
    Base object for Caches
    """
    def _hash(self, key):
        return str(key)

    def __contains__(self, key):
        key_hash = self._hash(key)
        return self._in(key_hash)

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
    def _in(self, key):
        raise NotImplementedError

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
        if self._in(key_hash):
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


class DictionaryCache(BaseCache):
    def __init__(self) -> None:
        super().__init__()
        self._cache = {}

    def _set(self, key, value):
        self._cache[key] = value

    def _get(self, key):
        return self._cache[key]

    def _del(self, key):
        del self._cache[key]

    def _in(self, key):
        return key in self._cache


class MemoryCache(DictionaryCache):
    """
    A dictionary in-memory cache
    """
    def load(self):
        pass

    def save(self):
        pass


class JsonCache(DictionaryCache):
    """
    A cache for downloaded TLE data using a json file
    """
    filename = 'passpredict.json'

    def __init__(self, filename=filename):
        super().__init__()
        self.filename = filename

    def load(self, strict=False):
        """ Load cache from json file """
        try:
            with open(self.filename, 'r') as f:
                self._cache = json.load(f)
        except FileNotFoundError as e:
            if strict:
                raise e
            else:
                pass

    def save(self):
        """ Save cache to json file """
        with open(self.filename, 'w') as f:
            json.dump(self._cache, f, indent=2)


class ShelfCache(DictionaryCache):
    """
    A cache for downloaded TLE data using a shelf.
    """
    filename = 'passpredict.db'

    def __init__(self, filename=filename):
        super().__init__()
        self.filename = filename

    def load(self, strict=False):
        """  Load shelf db from file  """
        if strict:
            flag = 'r'
        else:
            # Create the shelf if it doesn't exist
            flag = 'c'
        self._cache = shelve.DbfilenameShelf(self.filename, flag=flag)
        return self

    def save(self):
        self._cache.close()
