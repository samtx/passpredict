# TLE source for web app
import abc
import datetime
from typing import NamedTuple, Tuple

from orbit_predictor.sources import MemoryTLESource

from .satellites import SatellitePredictor
from .base import TLESource

# from orbit_predictor.sources
class TLE(NamedTuple):
    sate_id: int        # NORAD satellite ID
    lines: Tuple[str]   # tuple of tle strings (tle1, tle2)
    date: datetime.datetime   # datetime in UTC


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
