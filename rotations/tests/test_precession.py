# test_precession.py

import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose
import pytest

from .. import precession 
from ..core import mxv, mxmxm
import passpredict.constants as constants
from passpredict.utils import epoch_from_tle_datetime
from passpredict.timefn import julian_date, jd2jc
from passpredict.constants import DEG2RAD, RAD2DEG


@pytest.mark.parametrize(
    "epoch_string, zeta_expected, theta_expected, z_expected",
    [
        ('00182.78495062', 0.0031796*DEG2RAD, 0.0027633*DEG2RAD, 0.0031796*DEG2RAD),  # Vallado, p. 234
        ('00179.78495062', 0.0031270*DEG2RAD, 0.0027176*DEG2RAD, 0.0031270*DEG2RAD),  # Vallado, p. 234
    ],
    ids=[
        'epoch=00182.78495062',
        'epoch=00179.78495062',
    ]
)
def test_fk5_precession_angles_from_epoch(epoch_string, zeta_expected, theta_expected, z_expected):
    """
    Vallado, p.234
    """
    epoch = epoch_from_tle_datetime(epoch_string)
    jd = julian_date(epoch)
    tt = jd2jc(jd)
    zeta, theta, z = precession.fk5_precession_angles(tt)
    assert_almost_equal(zeta, zeta_expected, decimal=9)
    assert_almost_equal(theta, theta_expected, decimal=9)
    assert_almost_equal(z, z_expected, decimal=9)


def test_fk5_precession_angles():
    """
    Vallado, Eg. 3-15, p.231
    """
    # April 6, 2004, 07:51:28.386009 UTC
    tt = 0.0426236319  # Julian centuries since J2000
    zeta, theta, z = precession.fk5_precession_angles(tt)
    assert_almost_equal(zeta, 0.0273055*DEG2RAD, decimal=9)
    assert_almost_equal(theta, 0.0237306*DEG2RAD, decimal=9)
    assert_almost_equal(z, 0.0273059*DEG2RAD, decimal=9)


def test_precession_rotation():
    """
    Vallado, p.231
    
    Rotate from MOD --> GCRF
    """
    rTOD = np.array((5094.0283745, 6127.8708164, 6380.2485164))
    rGCRF_expected = np.array((5102.508958, 6123.011401, 6378.136928))
    tt = 0.0426236319
    prec = precession.fk5_precession(tt)
    rGCRF = mxv(prec.mtx, rTOD)
    assert_almost_equal(rGCRF[0], rGCRF_expected[0], decimal=6)
    assert_almost_equal(rGCRF[1], rGCRF_expected[1], decimal=6)
    assert_almost_equal(rGCRF[2], rGCRF_expected[2], decimal=6)


# def test_fk5_precession():
#     """
#     Vallado, Eg. 3-15, p.231
#     """
#     tt = 0.0426236319
#     rMOD = np.array([5094.0283745, 6127.8708164, 6380.2485164])
#     # zeta = 0.0273055 * constants.DEG2RAD
#     # theta = 0.0237306 * constants.DEG2RAD
#     # z = 0.0273059 * constants.DEG2RAD
#     prec = precession.fk5_precession(tt)
#     rGCRF_true = np.array([5102.508958, 6123.011401, 6378.136928])
#     assert_allclose(r, rGCRF_true)



if __name__ == "__main__":
    pytest.main(['-v', __file__])