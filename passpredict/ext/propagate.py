import numpy as np
from sgp4.api import Satrec, WGS84

from app.astrodynamics._rotations import teme2ecef
from app.astrodynamics._solar import sun_sat_illumination_distance
from app.schemas import Tle
from app.models import SatPredictData

# import debugpy
# debugpy.debug_this_thread()


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
    # r = np.reshape(r.ravel(order='F'), (3, r.shape[0]))
    # v = np.reshape(v.ravel(order='F'), (3, v.shape[0]))
    return r, v


def compute_satellite_data(tle: Tle, jd: np.ndarray, sun_rECEF: np.ndarray = None) -> SatPredictData:
    """
    Compute satellite data for Time

    Reference:
        https://docs.astropy.org/en/latest/coordinates/satellites.html

    """
    rTEME, _ = propagate_satellite(tle.tle1, tle.tle2, jd)
    rECEF = teme2ecef(jd, rTEME)
    if sun_rECEF is not None:
        sun_sat_dist = sun_sat_illumination_distance(rECEF, sun_rECEF)
        satpredictdata = SatPredictData(
            id=tle.satid, rECEF=rECEF, sun_sat_dist=sun_sat_dist
        )
    else:
        satpredictdata = SatPredictData(id=tle.satid, rECEF=rECEF)
    return satpredictdata
