# test rotations.py
import datetime

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

from passpredict import rotations
from passpredict import topocentric
from passpredict.models import RhoVector


def test_ecef2sez():
    """
    Vallado, Eg. 11-6, p.912
    """
    phi = 42.38  # latitude, deg
    lmda = -71.13  # longitude, deg
    # lmda = 136.2944
    h = 24  # height, m
    rsat = np.array([885.7296, -4389.3856, 5070.1765])
    rsite = topocentric.site_ECEF2(phi, lmda, h)
    rhoECEF = rsat - rsite
    print(rhoECEF)
    rSEZ = rotations.ecef2sez(rhoECEF, phi, lmda)
    rSEZ_true = np.array([-773.8654, -581.4980, 328.8145])
    np.set_printoptions(precision=8)
    assert_allclose(rSEZ, rSEZ_true)
    # for i in [0, 1, 2]:
    #     assert_almost_equal(rSEZ[i], rSEZ_true[i], decimal=0, verbose=True)


if __name__ == "__main__":
    import pytest
    pytest.main(['-v', __file__])