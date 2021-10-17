# TLE source for web app
from __future__ import annotations
import abc
import datetime
from typing import NamedTuple, Tuple

from orbit_predictor.sources import (
    TLESource,
    MemoryTLESource as MemoryTLESourceBase,
)

from typing import TYPE_CHECKING

from passpredict.predictors import SatellitePredictor
if TYPE_CHECKING:
    from .predictors import SatellitePredictorBase


# from orbit_predictor.sources
class TLE(NamedTuple):
    sate_id: int        # NORAD satellite ID
    lines: Tuple[str]   # tuple of tle strings (tle1, tle2)
    date: datetime.datetime   # datetime in UTC


class MemoryTLESource(MemoryTLESourceBase):

    def get_predictor(self, satid, predictor_class: SatellitePredictorBase = SatellitePredictor) -> SatellitePredictorBase:
        klass = predictor_class(satid, self)
        return klass


class AsyncTLESourceBase(TLESource):
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
    async def get_predictor(self, satid: int, date: datetime.datetime) -> SatellitePredictorBase:
        """
        Create Predictor instance with TLE data
        """
        raise NotImplementedError



