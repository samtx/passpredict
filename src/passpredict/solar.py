from functools import lru_cache
import warnings

from . import _solar
from .constants import MJD0


@lru_cache(maxsize=128)
def sun_pos(jd: float):
    warnings.warn("sun_pos() will take MJD times in future releases", DeprecationWarning)
    mjd = jd - MJD0
    return _solar.sun_pos_mjd(mjd)


@lru_cache(maxsize=128)
def sun_pos_mjd(mjd: float):
    return _solar.sun_pos_mjd(mjd)
