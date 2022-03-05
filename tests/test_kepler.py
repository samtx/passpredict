from math import radians, sqrt, degrees
import datetime

import numpy as np
import pytest
from pytest import approx

from passpredict.orbit import kepler_eqn_E, nu_to_anomaly, anomaly_to_nu, kepler_eqn
from passpredict.satellites.kepler import pkepler, pkepler2, KeplerPropagator
from passpredict.time import julian_date_from_datetime


@pytest.mark.parametrize(
    'M, e, E_expected',
    [
        pytest.param(235.4, 0.4, 3.84866174509717, id='Vallado, Eg 2-1, p.66'),
    ]
)
def test_kepler_eqn_E(M, e, E_expected):
    Mrad = radians(M)
    E = kepler_eqn_E(Mrad, e)
    assert E == approx(E_expected, abs=1e-8)


@pytest.mark.parametrize(
    'e, E, nu',
    # test cases from orbit_predictor/tests/test_angles.py
    [
        # (0.0, 0.0, 0.0),
        (0.05, 10.52321, 11.05994),
        (0.10, 54.67466, 59.49810),
        (0.35, 142.27123, 153.32411),
        (0.61, 161.87359, 171.02189),
    ]
)
class TestTrueEccentricAnomalies:

    def test_nu_to_anomaly(self, nu, e, E):
        Ecomputed = nu_to_anomaly(radians(nu), e)
        Edeg = degrees(Ecomputed)
        assert Edeg == approx(E)
        nucomputed = anomaly_to_nu(Ecomputed, e)
        nudeg = degrees(nucomputed)
        assert nudeg == approx(nu)

    def test_anomaly_to_nu(self, nu, e, E):
        nucomputed = anomaly_to_nu(radians(E), e)
        nudeg = degrees(nucomputed)
        assert nudeg == approx(nu)
        Ecomputed = nu_to_anomaly(nucomputed, e)
        Edeg = degrees(Ecomputed)
        assert Edeg == approx(E)


@pytest.mark.parametrize(
    'e, M, nu',
    # test cases from orbit_predictor/tests/test_angles.py
    [
        # (0.0, 0.0, 0.0),
        (0.05, 10.0, 11.06),
        (0.06, 30.0, 33.67),
        (0.04, 120.0, 123.87),
        (0.14, 65.0, 80.50),
        (0.19, 21.0, 30.94),
        (0.35, 65.0, 105.71),
        (0.48, 180.0, 180.0),
        (0.75, 125.0, 167.57)
    ]
)
class TestTrueMeanAnomalies:

    def test_nu_to_mean_anomaly(self, nu, e, M):
        E = nu_to_anomaly(radians(nu), e)
        Mcomputed = kepler_eqn(E, e)
        Mdeg = degrees(Mcomputed)
        assert Mdeg == approx(M, abs=1e-2)

    def test_mean_to_true_anomaly(self, nu, e, M):
        E = kepler_eqn_E(radians(M), e)
        nucomputed = anomaly_to_nu(E, e)
        nudeg = degrees(nucomputed)
        assert nudeg == approx(nu, abs=1e-2)

def test_pkepler():
    """
    Ref: orbit-predictor/tests/test_numerical_predictor.py
    """
    a = 6780
    ecc = 0.001
    inc = 28.5
    raan = 67.0
    argp = 355.0
    nu = 250.0
    dt = datetime.timedelta(hours=3).total_seconds()  # plus 3 hours from epoch
    r, v = pkepler(
        dt,
        radians(argp),
        ecc,
        radians(inc),
        radians(raan),
        a,
        radians(nu),
    )
    assert r == approx([2085.9287615146, -6009.5713894563, -2357.3802307070], rel=1e-2)
    assert v == approx([6.4787522759177, 3.2366136616580, -2.5063420188165], rel=1e-2)


def test_pkepler_from_rv():
    """
    Vallado, Eg 11-6, p.913
    """
    epoch_jd = 2450540.4
    r0 = np.array([6585.038266, 1568.184321, 9.116355])
    v0 = np.array([-1.1157766, 4.6316816, 6.0149576])
    ndot = 7.889e-5 * 2
    # ndot = 0
    satellite = KeplerPropagator.from_rv(r0, v0, epoch_jd, ndot=ndot)
    d = datetime.datetime(1997, 4, 2, 1, 8, tzinfo=datetime.timezone.utc)
    jd = sum(julian_date_from_datetime(d))
    dt = (jd - epoch_jd) * 86400.0  # dt in seconds
    argp0, e0, inc, raan0, a0, nu0, ndot, nddot = satellite.coe_init
    r, _ = pkepler(dt, argp0, e0, inc, raan0, a0, nu0, ndot=ndot, nddot=nddot)
    assert r == approx([-2811.2769, 3489.2632, 5069.5763])