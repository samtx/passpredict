# TLE source for web app
import abc
import datetime
import copy
from typing import NamedTuple, Tuple, Any, Union
from pathlib import Path
from collections import defaultdict

from orbit_predictor.sources import MemoryTLESource
import httpx

from .satellites import SatellitePredictor
from .base import TLESource
from .caches import JsonCache
from .tle import epoch_from_tle

# from orbit_predictor.sources
class TLE(NamedTuple):
    satid: Union[int, str]        # NORAD satellite ID
    lines: Tuple[str]   # tuple of tle strings (tle1, tle2)

    @property
    def sate_id(self):
        return self.satid

    @property
    def epoch(self) -> datetime:
        return epoch_from_tle(self.lines[0])


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
    def __init__(self) -> None:
        self.cache_path = Path.home() / '.passpredict' / 'tle.dat'
        self.data = {}

    def load(self):
        """
        Load TLEs from cache
        """
        with JsonCache(self.cache_path) as cache:
            self.data = copy.deepcopy(cache.cache)

    def save(self):
        """
        Save TLEs into cache file
        """
        with JsonCache(self.cache_path) as cache:
            for satkey in self.data:
                if "sat:" in satkey:
                    cache.set(f"sat:")

    def add_tle(self, satid: int, tle: TLE, epoch: datetime.datetime):
        """
        Add TLE to local cache
        """
        return super().add_tle(satid, tle, epoch)

    def _get_tle(self, satid: int, date):
        """
        Query TLE from Cache
        """
        return super()._get_tle(satid, date)

    def get_tle(self, satid: int, date):
        """
        Query TLE
        """
        return super().get_tle(satid, date)

    def get_predictor(self, satid: int):
        """
        Return satellite predictor from tle
        """
        predictor = SatellitePredictor(satid, self)
        predictor.set_propagator()
        return predictor


class MemoryTLESource(PasspredictTLESource):
    def __init__(self):
        self.tles = defaultdict()

    def add_tle(self, tle: TLE):
        self.tles[tle.satid] = tle

    def get_tle(self, satid: Union[int, str]):
        return self.tles[satid]


