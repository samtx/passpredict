import datetime

import pytest

from passpredict import Location, Observer, SGP4Propagator, TLE


@pytest.fixture()
def observer():
    satid = 'ISS'
    tle_lines = (
        "1 25544U 98067A   20196.51422950 -.00000046  00000-0  72206-5 0  9999",
        "2 25544  51.6443 213.2207 0001423 114.8006 342.8278 15.49514729236251"
    )
    tle = TLE(satid, tle_lines)
    satellite = SGP4Propagator.from_tle(tle)
    location = Location("Austin, Texas", 30.2711, -97.7434, 0)
    return Observer(location, satellite)


@pytest.mark.parametrize('method', ('op', 'brute',))
@pytest.mark.parametrize('tol', (0.5, 0.25, 0.1,))
@pytest.mark.parametrize('time_step', (60, 20, 5,))
def test_observer_benchmark_method(observer, method, tol, time_step):
    if method != "brute" and time_step != 60:
        # we don't need to run unit test on different
        # time steps for non brute force methods
        return
    start = datetime.datetime(2020, 7, 14, tzinfo=datetime.timezone.utc)
    end = start + datetime.timedelta(days=10)
    overpasses = observer.pass_list(
        start_date=start,
        limit_date=end,
        method=method,
        aos_at_dg=10,
        tol=tol,
    )
    assert len(overpasses) > 0

