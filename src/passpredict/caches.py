import json
import shelve
import time
import pathlib

class JsonCache:
    """
    A cache for downloaded TLE data using a json file
    """
    filename = 'tle.json'

    def __init__(self, filename=filename):
        # if not pathlib.Path(filename).is_file():
        #     raise Exception(f"JsonCache filename {filename} is not valid.")
        self.filename = filename
        self.cache = {}

    def __contains__(self, key):
        return key in self.cache

    def __setitem__(self, key, value):
        self.cache[key] = value

    def __getitem__(self, key):
        return self.cache[key]

    def __delitem__(self, key):
        del self.cache[key]

    def __enter__(self):
        self.load()
        return self

    def __exit__(self, *args):
        self.save()

    def load(self):
        """ Load cache from json file """
        try:
            with open(self.filename, 'r') as f:
                self.cache = json.load(f)
        except FileNotFoundError:
            pass

    def save(self):
        """ Save cache to json file """
        with open(self.filename, 'w') as f:
            json.dump(self.cache, f, indent=2)

    def get(self, key, ignore_ttl=False):
        if key not in self: return None
        item = self.cache[key]
        if (not ignore_ttl) and (ttl := item['ttl']):
            if time.time() > ttl:
                del self.cache[key]
                return None
        return item['data']

    def set(self, key, value, ttl=None):
        if ttl is not None:
            ttl += time.time()
        item = {'ttl': ttl, 'data': value}
        self[key] = item


class ShelfCache:
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
