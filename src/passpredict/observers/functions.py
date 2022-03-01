from __future__ import annotations
import typing
import datetime

import numpy as np

from ..time import julian_date


def make_utc(d: datetime.datetime) -> datetime.datetime:
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


def find_min(f: typing.Callable, a: float, b: float, tol: float) -> float:
    """
    Find minimum of bounded univariate scalar function using scipy.optimize.minimize_scalar
    """
    assert a < b
    # res = minimize_scalar(f, bounds=(a, b), method='golden', tol=tol, options={'xatol': tol})
    diff = b - a
    fvec = np.vectorize(f)
    N = 5
    while diff > tol:
        x = np.linspace(a, b, N)
        y = fvec(x)
        i = y.argmin()
        if i == N:
            a = x[N - 1]
            b = x[N]
        elif i == 0:
            a = x[0]
            b = x[1]
        else:
            a = x[i - 1]
            b = x[i + 1]
        diff = abs(b - a)
    xsol = a + diff / 2
    fsol = f(xsol)
    return xsol, fsol


def julian_date_sum(d: datetime.datetime) -> float:
    jd, jdfr = julian_date(d.year, d.month, d.day, d.hour, d.minute, d.second + d.microsecond/1e6)
    return jd + jdfr
