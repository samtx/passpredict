from itertools import zip_longest
import datetime
from math import floor

from numpy import pi as np_pi


def shift_angle(x: float) -> float:
    """Shift angle in radians to [-pi, pi)

    Args:
        x: float, angle in radians

    Reference:
        https://stackoverflow.com/questions/15927755/opposite-of-numpy-unwrap/32266181#32266181
    """
    return (x + np_pi) % (2 * np_pi) - np_pi


def grouper(iterable, n, fillvalue=None):
    """
    from itertools recipes https://docs.python.org/3.7/library/itertools.html#itertools-recipes
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def datetime_now_utc():
    """
    Get datetime.datetime object in current utc
    """
    return datetime.datetime.now(datetime.timezone.utc)


def direction_from_azimuth(azimuth: float) -> str:
    ''' Return ordinal direction from azimuth degree '''
    azm = azimuth % 360
    mod = 360/16. # number of degrees per coordinate heading
    start = 0 - mod/2
    n = int(floor((azm-start)/mod))
    coordinates = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','N']
    return coordinates[n]