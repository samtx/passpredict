import numpy as np
from astropy import units as u
from astropy.coordinates import TEME, CartesianRepresentation, ITRS
from astropy.time import Time
from sgp4.api import Satrec, WGS84

from passpredict.solar import is_sat_illuminated
from passpredict.schemas import Tle
from passpredict.models import Sat, SatPredictData, SunPredictData


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


def compute_satellite_data(tle: Tle, t: Time, sun: SunPredictData = None) -> SatPredictData:
    """
    Compute satellite data for Time
    
    Reference: 
        https://docs.astropy.org/en/latest/coordinates/satellites.html

    """
    sat = Sat()
    sat.time = t
    r, _ = propagate_satellite(tle.tle1, tle.tle2, t.jd)
    # Use the TEME reference frame from astropy
    teme = TEME(CartesianRepresentation(r * u.km), obstime=t)
    ecef = teme.transform_to(ITRS(obstime=t))
    sat.rECEF = ecef.data.xyz.value.astype(np.float32)  # extract numpy array from astropy object
    # sat.subpoint = ecef.earth_location
    # sat.latitude = sat.subpoint.lat.value
    # sat.longitude = sat.subpoint.lon.value
    if sun is not None:
        sat.illuminated = is_sat_illuminated(sat.rECEF, sun.rECEF)
    return SatPredictData(id=tle.satid, rECEF=sat.rECEF, illuminated=sat.illuminated)