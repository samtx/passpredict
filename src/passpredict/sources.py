# TLE source for web app
from __future__ import annotations
import abc
import datetime
import copy
from typing import Union, TYPE_CHECKING, List
from pathlib import Path
from collections import defaultdict

import httpx

from .satellites import SGP4Propagator
from .caches import BaseCache, MemoryCache
from .exceptions import CelestrakError
from .orbit import TLE
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
        predictor = SGP4Propagator(satid, self)
        predictor.tle = tle
        predictor.set_propagator()
        return predictor


class CelestrakTLESource(PasspredictTLESource):
    """
    TLE source that checks the cache for orbital elements,
    otherwise queries Celestrak website
    """
    def __init__(self, cache: BaseCache = None) -> None:
        if cache:
            self.cache = cache
        else:
            self.cache = MemoryCache()

    def add_tle(self, satid: int, tle: TLE, epoch: datetime.datetime):
        """
        Add TLE to local cache
        """
        return super().add_tle(satid, tle, epoch)

    def get_tle(self, satid: int) -> TLE:
        """
        Query TLE from cache. If not found in cache, then query Celestrak url
        """
        key = f"tle:{satid}"
        res = self.cache.get(key)
        if res:
            tle = TLE(res['satid'], res['lines'], name=res.get('name', ''))
        else:
            tle = self._query_tle_from_celestrak(satid)
            self.cache.set(key, tle.dict(), ttl=86400)
        return tle

    def get_tle_category(self, category: str) -> List[TLE]:
        """
        Query TLE from cache. If not found in cache, then query Celestrak url
        """
        # first, get list of satellite IDs for category
        key = f"tlecategory:{category}"
        res = self.cache.get(key)
        if res:
            satids = res
            tles = []
            for satid in satids:
                tle = self.get_tle(satid)
                tles.append(tle)
        else:
            tles = self._query_tle_category_from_celestrak(category)
            satid_list = [t.satid for t in tles]
            self.cache.set(key, satid_list, ttl=86400)
            for tle in tles:
                self.cache.set(f"tle:{tle.satid}", tle.dict(), ttl=86460)
        return tles

    def _query_tle_from_celestrak(self, satid: int = None):
        """
        Download current TLEs from Celestrak and save them to a JSON file
        """
        url = 'https://celestrak.com/NORAD/elements/gp.php'
        params = {
            'CATNR': satid,
            'FORMAT': 'TLE',
        }
        r = httpx.get(url, params=params)
        if r.text.lower() in ("no tle found", "no gp data found") or r.status_code >= 300:
            raise CelestrakError(f'Celestrak TLE for satellite {satid} not found')
        tle_strings = r.text.splitlines()
        tle = parse_tle(tle_strings)
        return tle

    def _query_tle_category_from_celestrak(self, category: str) -> List[TLE]:
        """
        Download current TLEs from Celestrak and save them to a JSON file

        Enter Celestrak category string.
        Eg. visual, stations, tle-new, weather, noaa, oneweb, starlink
        """
        url = f'https://celestrak.com/NORAD/elements/{category}.txt'
        r = httpx.get(url)
        if not r.text or r.status_code >= 300:
            raise CelestrakError(f'Celestrak TLEs for {category}.txt not found')
        tle_strings = r.text.splitlines()
        tles = parse_multiple_tles(tle_strings)
        return tles


class MemoryTLESource(PasspredictTLESource):
    def __init__(self):
        self.tles = defaultdict()

    def add_tle(self, tle: TLE):
        self.tles[tle.satid] = tle

    def get_tle(self, satid: Union[int, str]):
        return self.tles[satid]


def parse_multiple_tles(tle_lines):
    """ Assumes TLEs are 3 line TLEs """
    tles = []
    for raw_tle in grouper(tle_lines, 3):
        tle = parse_tle(raw_tle)
        tles.append(tle)
    return tles


def parse_tle(tle_lines):
    """
    Parse a single 2-line or 3-line TLE from celestrak
    """
    if len(tle_lines) == 2:
        tle1, tle2 = tle_lines
        name = ""
    elif len(tle_lines) == 3:
        tle0, tle1, tle2 = tle_lines
        name = tle0.strip()  # satellite name
    else:
        raise Exception(f"Invalid TLE strings {tle_lines}")
    satid = int(tle1[2:7])
    tle = TLE(satid, (tle1, tle2), name=name)
    return tle