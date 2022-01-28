# TLE source for web app
from __future__ import annotations
import abc
import datetime
import copy
from typing import Union, TYPE_CHECKING
from pathlib import Path
from collections import defaultdict

import httpx

from .satellites import SatellitePredictor
from .base import TLESource
from .caches import JsonCache
from .tle import TLE
from .utils import grouper


class PasspredictTLESource(abc.ABC):

    @abc.abstractmethod
    def get_tle(self, satid: int):
        raise NotImplementedError

    @abc.abstractmethod
    def add_tle(self, satid: int, tle: TLE):
        """
        Add TLE to source object
        """
        raise NotImplementedError


class AsyncPasspredictTLESource(abc.ABC):
    """
    TLE source that checks the redis cache and postgres database
    for orbital elements
    """

    @abc.abstractmethod
    async def add_tle(self, satid, tle):
        """
        Add TLE to cache
        """
        raise NotImplementedError

    @abc.abstractmethod
    async def get_tle_or_404(self, satid: int, date: datetime.datetime):
        """
        Get TLE from source object
        Raise 404 error if TLE not found in database
        """
        raise NotImplementedError

    @abc.abstractmethod
    async def get_predictor(self, satid: int, date: datetime.datetime):
        """
        Create Predictor instance with TLE data
        """
        tle = await self.get_tle_or_404(satid, date)
        predictor = SatellitePredictor(satid, self)
        predictor.tle = tle
        predictor.set_propagator()
        return predictor


class CelestrakTLESource(TLESource):
    """
    TLE source that checks the cache for orbital elements,
    otherwise queries Celestrak website
    """
    def __init__(self, cache_path: Union[str, Path] = None) -> None:
        if not cache_path:
            cache_path = 'tle.json'
        self.cache_path = cache_path
        self.cache = JsonCache(self.cache_path)

    def load(self, strict: bool = False):
        """
        Load TLEs from cache
        """
        try:
            self.cache.load()
        except FileNotFoundError as e:
            if strict:
                raise e
            else:
                pass

    def save(self):
        """
        Save TLEs into cache file
        """
        self.cache.save()

    def add_tle(self, satid: int, tle: TLE, epoch: datetime.datetime):
        """
        Add TLE to local cache
        """
        return super().add_tle(satid, tle, epoch)

    def get_tle(self, satid: int):
        """
        Query TLE from cache. If not found in cache, then query Celestrak url
        """
        key = f"tle:{satid}"
        res = self.cache.get(key)
        if res:
            tle = TLE(res['satid'], res['lines'], name=res.get('name', ''))
        else:
            tle = self._query_tles_from_celestrak(satid)
            self.cache.set(key, tle.dict(), ttl=86400)
        return tle

    def _query_tles_from_celestrak(self, satid: int = None):
        """
        Download current TLEs from Celestrak and save them to a JSON file
        """
        url = 'https://celestrak.com/satcat/tle.php'
        params = {'CATNR': satid}
        r = httpx.get(url, params=params)
        tle_strings = r.text.splitlines()
        tle = self._parse_tle(tle_strings)
        return tle

    def _parse_tle(self, tle_strings):
        """
        Parse a single 3-line TLE from celestrak
        """
        if len(tle_strings) == 2:
            tle1, tle2 = tle_strings
            name = ""
        elif len(tle_strings) == 3:
            tle0, tle1, tle2 = tle_strings
            name = tle0.strip()  # satellite name
        else:
            raise Exception(f"Invalid TLE strings {tle_strings}")
        satid = int(tle1[2:7])
        tle = TLE(satid, (tle1, tle2), name=name)
        return tle


class MemoryTLESource(PasspredictTLESource):
    def __init__(self):
        self.tles = defaultdict()

    def add_tle(self, tle: TLE):
        self.tles[tle.satid] = tle

    def get_tle(self, satid: Union[int, str]):
        return self.tles[satid]
