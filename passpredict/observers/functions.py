from __future__ import annotations
import typing
import datetime

from ..time import julian_date


def _make_utc(d: datetime.datetime) -> datetime.datetime:
    """ Make a datetime a UTC timezone aware datetime """
    if not d:
        return d
    if d.tzinfo:
        d = d.astimezone(datetime.timezone.utc)
    else:
        d = d.replace(tzinfo=datetime.timezone.utc)
    return d


def find_root(f: typing.Callable, a: float, b: float, tol: float) -> float:
    """
    Find root of function f() with bisection method within tolerance tol.
    Return root.
    """
    assert a < b
    fa, fb = f(a), f(b)
    if fa*fb >= 0:
        return None
        # # return value closest to zero
        # if abs(fa) < abs(fb):
        #     return a
        # else:
        #     return b

    diff = tol + 1
    while diff > tol:
        mid = (a + b) / 2
        fmid = f(mid)
        if fa*fmid < 0:
            b = mid
            fb = f(b)
        elif fb*fmid < 0:
            a = mid
            fa = f(a)
        diff = abs(b - a)
    return mid


def julian_date_sum(d: datetime.datetime) -> float:
    jd, jdfr = julian_date(d.year, d.month, d.day, d.hour, d.minute, d.second + d.microsecond/1e6)
    return jd + jdfr