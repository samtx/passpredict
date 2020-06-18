# test_solar.py
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import passpredict.timefn as timefn
import passpredict.solar as solar
from passpredict.constants import DEG2RAD, RAD2DEG, AU_KM
from datetime import datetime, timedelta


def test_sun_pos():
    """
    Vallado, Eg. 5-1, p. 280
    """
    dt = datetime(2006, 4, 2)  # April 2, 2006, 00:00 UTC
    jdt = timefn.julian_date2(dt)
    assert_almost_equal(jdt, 2453827.5, decimal=12)
    jdt = np.asarray(jdt)
    r = solar.sun_pos(jdt)
    r_true = np.array([146186212, 28788976, 12481064], dtype=np.float)
    r_true = np.reshape(r_true, (3, 1))
    assert_allclose(r, r_true, rtol=1e-4)


def test_sun_pos_2():
    """
    Vallado, Eg. 11-6, p. 913
    """
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = timefn.julian_date2(dt)
    jdt = np.asarray(jdt)
    r = solar.sun_pos(jdt) / AU_KM
    r_true = np.array([0.9765, 0.1960, 0.0850], dtype=np.float)
    r_true = np.reshape(r_true, (3, 1))
    assert_allclose(r, r_true, rtol=1e-3)


def test_sun_sat_angle():
    """
    Vallado, Eg. 11-6, p.913
    """
    dt = datetime(1997, 4, 2, 1, 8)  # April 2, 1997, 01:08:0.00 UTC
    jdt = timefn.julian_date2(dt)
    jdt = np.asarray(jdt)
    rsun = solar.sun_pos(jdt)
    rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
    rsun = np.atleast_2d(rsun)
    rsat = np.atleast_2d(rsat).T
    sunangle = solar.sun_sat_angle(rsat, rsun) * RAD2DEG
    assert_almost_equal(sunangle, 76.0407, decimal=3)


def test_sun_sat_angle2():
    """
    Vallado, Eg. 11-6, p.913
    """
    rsat = np.array([-2811.2769, 3486.2632, 5069.5763])
    rsun = np.array([0.9765, 0.1960, 0.0850]) * AU_KM
    sunangle = solar.sun_sat_angle(rsat, rsun) * RAD2DEG
    assert_almost_equal(sunangle, 76.0407, decimal=3)


def test_sun_sat_orthogonal_distance():
    """
    Vallado, Eg. 11-6, p.913
    """
    r = np.array([-2811.2769, 3486.2632, 5069.5763])  # sat, ECI coordinates
    zeta = 76.0407  # deg
    dist = solar.sun_sat_orthogonal_distance(r, zeta * DEG2RAD)
    assert_almost_equal(dist, 6564.6870, decimal=4)


