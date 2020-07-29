import datetime

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

from passpredict.rotations import transform
from passpredict.rotations import sidereal
from passpredict.rotations import nutation
from passpredict import constants
from passpredict import timefn
from passpredict import topocentric
from passpredict.utils import epoch_from_tle_datetime
from passpredict.constants import DEG2RAD, RAD2DEG


def test_thetaGMST_from_skyfield():
    """From skyfield.tests.test_earth_satellites.py"""
    ms = int(1e6) + 386000 - 439961
    dt = datetime.datetime(2004, 4, 6, 7, 51, 27, ms)
    jd = timefn.julian_date(dt)
    theta, thetadt = sidereal.theta_GMST1982(jd)
    assert_almost_equal(theta, 5.459562584754709, decimal=15)


@pytest.mark.parametrize(
    "dpsi, meaneps, eq_expected",
    [
        (-0.004250260*DEG2RAD, 23.43922657*DEG2RAD, -0.0038995*DEG2RAD),
        (-0.004337544*DEG2RAD, 23.43922657*DEG2RAD, -0.0039796*DEG2RAD),
    ],
    ids=[
        'epoch=00182.78495062',
        'epoch=00179.78495062',
    ]
)
def test_equinox1982_geometric_terms(dpsi, meaneps, eq_expected):
    """
    Vallado, p.234, epoch datetime = 00179.78495062
    """
    eq = sidereal.equinox1982_geometric_terms(dpsi, meaneps)
    assert_almost_equal(eq, eq_expected)


@pytest.mark.xfail(reason='Need to fix this later')
def test_equinox1982_geometric_terms_3():
    """
    Vallado, Eg.3-15, p.230,
    """
    theta_GMST1982 = 312.8098943  # degrees
    eq_true =  312.8067654 - theta_GMST1982
    jd = timefn.julian_date(2006, 4, 6, 7, 51, 28.386009)
    tt = timefn.jd2jc(jd)
    nut = nutation.fk5_nutation(tt)
    eq = sidereal.equinox1982(nut.dpsi, nut.eps, nut.omega_moon)
    assert_almost_equal(eq, eq_true * constants.DEG2RAD)


@pytest.mark.skip('failing')
def test_theta_GAST1982():
    """
    Vallado, Eg. 3-15, p. 230
    """
    dut1 = -0.4399619
    jdut1 = timefn.julian_date(2004, 4, 6, 7, 51, 28.386009 + dut1)
    ttut1 = timefn.jd2jc(jdut1)
    gmst = 312.8098943 * constants.DEG2RAD
    tt = 0.0426236319
    nut = nutation.fk5_nutation(tt)
    print('dpsi ',nut.dpsi*constants.RAD2DEG)
    print('eps ',nut.eps*constants.RAD2DEG)
    print('om ',nut.omega_moon % 360.0)
    # eps = 23.4387368 * constants.DEG2RAD
    eps = 23.4407685 * constants.DEG2RAD
    dpsi = -0.0034108 * constants.DEG2RAD
    omega_moon = 42.6046140 * constants.DEG2RAD
    Eq = sidereal.equinox1982(dpsi, eps, omega_moon)
    gast = sidereal.theta_GAST1982(gmst, Eq)
    assert_almost_equal(gast * constants.RAD2DEG, 312.8067654)



def test_GMST():
    """
    Vallado, p.230
    """
    tt_ut1 = 0.0426236319  # julian centuries for UT1
    gmst_expected = 312.8098943  # degrees
    gmst = sidereal.theta_GMST1982(tt_ut1)

def test_sidereal_rotation():
    """
    Vallado, p.230

    PEF -> TOD rotation
    """
    tt = 0.0426236319
    rPEF = np.array((-1033.4750313, 7901.3055856, 6380.3445328))
    rTOD_expected = np.array((5094.5147804, 6127.3664612, 6380.3445328))

