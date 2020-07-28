import datetime
from pathlib import Path
import os

import numpy as np

from passpredict import propagate
from passpredict.rotations.transform import teme2ecef
from passpredict.rotations.polar import eop
from passpredict.timefn import julian_date
from passpredict.predictions import compute_time_array

tz_utc = datetime.timezone.utc

def test_propagate_iss():
    """
    Compare results of propagate() with skyfield EarthSatellite.at()

    datetime_start = 2020-06-01:00:00:00 UTC
    datetime_end   = 2020-06-11:00:00:00 UTC
    dt_seconds = 5
    
    data/skyfield_iss_rECEF.npy
    """
    tle1 = "1 25544U 98067A   20154.57277630  .00016717  00000-0  10270-3 0  9118"
    tle2 = "2 25544  51.6443  60.8122 0001995  12.6931 347.4269 15.49438452 29742"
    datetime_start = datetime.datetime(2020, 6, 1, 0, 0, 0, tzinfo=tz_utc)
    datetime_end = datetime.datetime(2020, 6, 11, 0, 0, 0, tzinfo=tz_utc)
    dt_sec = 5
    t = compute_time_array(datetime_start, datetime_end, dt_sec)
    rTEME, _ = propagate.propagate_satellite(tle1, tle2, jd=t.jd)#, dt0=datetime_start, dtf=datetime_end, dtsec=dt_sec)
    dUTC1, xp, yp = eop(t.jd)
    t.jd_utc1 = t.jd + dUTC1
    passpredict_rECEF = teme2ecef(rTEME, t.jd_utc1, xp, yp)    
    fname = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/skyfield_iss_rECEF.npy'))
    skyfield_rECEF = np.load(fname, allow_pickle=True)
    diff = np.linalg.norm(passpredict_rECEF - skyfield_rECEF, axis=0)

    # TODO: Make this more accurate in the future!!
    np.testing.assert_array_less(diff, 20.0, verbose=True)  # difference less than 20 km


if __name__ == "__main__":
    import pytest
    pytest.main(['-v', __file__])