from __future__ import annotations
import datetime
import typing
from math import degrees, floor
from functools import cached_property
from dataclasses import dataclass
from enum import Enum

if typing.TYPE_CHECKING:
    from ..locations import Location


class RangeAzEl(typing.NamedTuple):
    range: float  # km
    az: float     # deg
    el: float     # deg


@dataclass(frozen=True)
class PassPoint:
    dt: datetime                # datetime UTC
    range: float                # km
    azimuth: float              # deg
    elevation: float            # deg
    brightness: float = None    # magnitude brightness

    @cached_property
    def direction(self) -> str:
        ''' Return ordinal direction from azimuth degree '''
        azm = self.azimuth % 360
        mod = 360/16. # number of degrees per coordinate heading
        start = 0 - mod/2
        n = int(floor((azm-start)/mod))
        coordinates = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','N']
        return coordinates[n]


class PassType(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'


class BasicPassInfo:
    """
    Holds basic pass information:
        aos datetime
        los datetime
        tca datetime
        max elevation
        duration (property)
        valid (property)
    """
    def __init__(
        self,
        aos_dt: datetime.datetime,
        tca_dt: datetime.datetime,
        los_dt: datetime.datetime,
        max_elevation: float,
        type_: PassType = None,
        vis_begin_dt: datetime.datetime = None,
        vis_end_dt: datetime.datetime = None,
    ):
        self.aos_dt = aos_dt if aos_dt is not None else None
        self.tca_dt = tca_dt
        self.los_dt = los_dt if los_dt is not None else None
        self.max_elevation = max_elevation
        self.type = type_
        self.vis_begin_dt = vis_begin_dt
        self.vis_end_dt = vis_end_dt

    @property
    def aos(self):
        """ Backwards compatibilty from orbit_predictor """
        return self.aos_dt

    @property
    def tca(self):
        """ Backwards compatibilty from orbit_predictor """
        return self.tca_dt

    @property
    def los(self):
        """ Backwards compatibilty from orbit_predictor """
        return self.los_dt

    @property
    def valid(self):
        if (self.aos is None) or (self.los is None):
            return False
        return (self.max_elevation > 0)

    @cached_property
    def max_elevation_deg(self):
        return degrees(self.max_elevation)

    @cached_property
    def duration(self) -> datetime.timedelta:
        return self.los - self.aos


@dataclass
class PredictedPass:
    satid: int
    location: Location
    aos: PassPoint
    tca: PassPoint
    los: PassPoint
    type: PassType = None
    azimuth: typing.Sequence[float] = None
    elevation: typing.Sequence[float] = None
    range: typing.Sequence[float] = None
    datetime: typing.Sequence[datetime.datetime] = None
    vis_begin: PassPoint = None
    vis_end: PassPoint = None
    brightness: float = None

    @cached_property
    def midpoint(self):
        """ Return datetime in UTC of midpoint of pass """
        midpt = self.aos.dt + (self.los.dt - self.aos.dt) / 2
        return midpt

    @cached_property
    def duration(self):
        """ Return pass duration in seconds """
        return (self.los.dt - self.aos.dt).total_seconds()

    def __repr__(self):
        return f"<PredictedPass {self.satid} over {repr(self.location)} on {self.aos.dt}"

