from itertools import zip_longest

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