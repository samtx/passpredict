from typing import NamedTuple
import datetime
from functools import lru_cache
from itertools import zip_longest

from timezonefinder import TimezoneFinder

try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo


tf = TimezoneFinder()


class DatetimeMetadata(NamedTuple):
    start: datetime.datetime
    end: datetime.datetime
    n_steps = int


def round_to_nearest_second(dt: datetime.datetime) -> datetime.datetime:
    """  Round datetime to nearest whole second  """
    if dt.microsecond >= 500_000:
        dt += datetime.timedelta(seconds=1)
    return dt.replace(microsecond=0)


def get_pass_detail_datetime_metadata(
    pass_,
    delta_s: float,
    *,
    pad_minutes: int = 5,
):
    # add cushion to either side, rounded to the nearest second
    padding = datetime.timedelta(minutes=pad_minutes)
    start_date = round_to_nearest_second(pass_.aos.dt - padding)
    end_date = round_to_nearest_second(pass_.los.dt + padding)
    n_steps_float = (end_date - start_date).total_seconds() / delta_s
    n_steps = round(n_steps_float)
    time_step = datetime.timedelta(seconds=delta_s)
    return start_date, n_steps, time_step


@lru_cache(maxsize=128)
def get_timezone_from_latlon(latitude: float, longitude: float) -> ZoneInfo:
    """
    Returns timezone ZoneInfo object
    """
    tz_str = tf.timezone_at(lng=longitude, lat=latitude)
    tz = datetime.timezone.utc if not tz_str else ZoneInfo(tz_str)
    return tz


def grouper(iterable, n, fillvalue=None):
    """
    from itertools recipes https://docs.python.org/3.7/library/itertools.html#itertools-recipes
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
