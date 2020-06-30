# test_precession.py

import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose
import pytest

import passpredict.precession as precession 
import passpredict.constants as constants



def test_fk5_precession_angles():
    """
    Vallado, Eg. 3-15, p.231
    """
    # April 6, 2004, 07:51:28.386009 UTC
    tt = 0.0426236319  # Julian centuries since J2000
    zeta, theta, z = precession.fk5_precession_angles(tt)
    assert_almost_equal(zeta, 0.0273055 * constants.DEG2RAD, decimal=9)
    assert_almost_equal(theta, 0.0237306 * constants.DEG2RAD, decimal=9)
    assert_almost_equal(z, 0.0273059 * constants.DEG2RAD, decimal=9)


def test_fk5_precession():
    """
    Vallado, Eg. 3-15, p.231
    """
    tt = 0.0426236319
    rMOD = np.array([5094.0283745, 6127.8708164, 6380.2485164])
    # zeta = 0.0273055 * constants.DEG2RAD
    # theta = 0.0237306 * constants.DEG2RAD
    # z = 0.0273059 * constants.DEG2RAD
    prec = precession.fk5_precession(tt)
    rGCRF_true = np.array([5102.508958, 6123.011401, 6378.136928])
    assert_allclose(rGCRF, rGCRF_true)



if __name__ == "__main__":
    pytest.main(['-v', __file__])