# TLE source for web app
import abc
import datetime
from typing import NamedTuple, Tuple
from pathlib import Path

from orbit_predictor.sources import MemoryTLESource
import httpx

from .satellites import SatellitePredictor
from .base import TLESource
from .caches import JsonCache

# from orbit_predictor.sources
class TLE(NamedTuple):
    satid: int        # NORAD satellite ID
    lines: Tuple[str]   # tuple of tle strings (tle1, tle2)
    date: datetime.datetime   # datetime in UTC

    @property
    def sate_id(self):
        return self.satid



class CelestrakTLESource(TLESource):
    """
    TLE source that checks the cache for orbital elements,
    otherwise queries Celestrak website
    """
    def __init__(self) -> None:
        self.cache_path = Path.home() / '.passpredict' / 'cache.dat'
        self.cache = JsonCache(self.cache_path)


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
        return super().get_predictor(satid)



class AsyncPasspredictTLESource(TLESource):
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
        predictor = SatellitePredictor(satid)
        predictor._source = self
        predictor.tle = tle
        predictor.set_propagator()
        return predictor
