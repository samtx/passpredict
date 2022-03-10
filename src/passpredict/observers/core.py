from __future__ import annotations
import datetime
from math import degrees, floor
from typing import TYPE_CHECKING, NamedTuple, Sequence
from functools import cached_property
from dataclasses import dataclass, asdict
from enum import Enum

if TYPE_CHECKING:
    from ..locations import Location


class RangeAzEl(NamedTuple):
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
    type: Visibility = None       # visibility status

    @cached_property
    def direction(self) -> str:
        ''' Return ordinal direction from azimuth degree '''
        azm = self.azimuth % 360
        mod = 360/16.0  # number of degrees per coordinate heading
        start = 0 - mod/2
        n = int(floor((azm-start)/mod))
        coordinates = [
            'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW',
            'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N'
        ]
        return coordinates[n]


class PassType(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'


class Visibility(str, Enum):
    daylight = 'daylight'
    unlit = 'unlit'
    visible = 'visible'
    not_visible = 'not visible'


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
        aos_mjd: float,
        tca_mjd: float,
        los_mjd: float,
        max_elevation: float,
        type_: PassType = None,
        vis_begin_mjd: float = None,
        vis_end_mjd: float = None,
        vis_tca_mjd: float = None,
    ):
        self.aos_mjd = aos_mjd
        self.tca_mjd = tca_mjd
        self.los_mjd = los_mjd
        self.max_elevation = max_elevation
        self.type = type_
        self.vis_begin_mjd = vis_begin_mjd
        self.vis_end_mjd = vis_end_mjd
        self.vis_tca_mjd = vis_tca_mjd

    @property
    def aos(self):
        return self.aos_mjd

    @property
    def tca(self):
        return self.tca_mjd

    @property
    def los(self):
        return self.los_mjd

    @cached_property
    def max_elevation_deg(self):
        return degrees(self.max_elevation)

    @cached_property
    def duration(self) -> float:
        """
        Returns duration of pass in seconds
        """
        return (self.los_mjd - self.aos_mjd) * 86400.0


@dataclass
class PredictedPass:
    satid: int
    location: Location
    aos: PassPoint
    tca: PassPoint
    los: PassPoint
    type: PassType = None
    azimuth: Sequence[float] = None
    elevation: Sequence[float] = None
    range: Sequence[float] = None
    datetime: Sequence[datetime.datetime] = None
    vis_begin: PassPoint = None
    vis_end: PassPoint = None
    vis_tca: PassPoint = None
    brightness: float = None

    @cached_property
    def midpoint(self):
        """ Return datetime in UTC of midpoint of pass """
        midpt = self.aos.dt + (self.los.dt - self.aos.dt) / 2
        return midpt

    @cached_property
    def duration(self) -> float:
        """ Return pass duration in seconds """
        return (self.los.dt - self.aos.dt).total_seconds()

    def __repr__(self):
        return (
            f"<PredictedPass {self.satid} over "
            f"{repr(self.location)} on {self.aos.dt}"
        )

    def dict(self):
        """
        Serialize into dictionary
        """
        data = {
            'satid': self.satid,
            'location': {
                'name': self.location.name,
                'lat': self.location.latitude_deg,
                'lon': self.location.longitude_deg,
                'h': self.location.elevation_m,
            },
            'aos': asdict(self.aos),
            'tca': asdict(self.tca),
            'los': asdict(self.los),
            'type': self.type
        }
        return data
