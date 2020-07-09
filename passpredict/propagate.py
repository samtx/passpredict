
import numpy as np
import functools
from passpredict.timefn import julian_date, jdt_tsince, invjday, jday2datetime, \
                                 jday2npdatetime64, truncate_datetime
from passpredict.solar import is_sat_illuminated, sun_pos
from passpredict.schemas import Tle
from passpredict.models import SatelliteRV
from passpredict.rotations.transform import teme2ecef
from passpredict.utils import parse_tles_from_celestrak, epoch_from_tle
import datetime
from sgp4.api import Satrec, WGS84
import json
import requests
import os

# use_cython = False
# try:
#     from passpredict._sgp4 import sgp4 as sgp4_pyx
#     use_cython = True
# except ImportError:
#     pass
# from passpredict.sgp4 import sgp4


def propagate_satellite(tle1, tle2, jd, *args, **kwargs):
    """Propagate satellite position forward in time.

    Parameters:
        tle1 : str
            first line of two line element set
        tle2 : str
            second line of two line element set
        jd : float (n)
            array of julian dates to propagate

    Returns:
        r : float (3, n)
            satellite position vector in TEME coordinates
        v : float (3, n)
            satellite velocity vector in TEME coordinates
    """

    satrec = Satrec.twoline2rv(tle1, tle2, WGS84)

    jd_array, fr_array = np.divmod(jd, 1)

    error, r, v = satrec.sgp4_array(jd_array, fr_array)

    # Change arrays from column major to row major while keeping C-continuous
    r = np.reshape(r.ravel(order='F'), (3, r.shape[0]))
    v = np.reshape(v.ravel(order='F'), (3, v.shape[0]))

    return r, v