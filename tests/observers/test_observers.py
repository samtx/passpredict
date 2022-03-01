from datetime import datetime, timedelta, timezone

import pytest

from passpredict import MemoryTLESource, Observer, Location
from passpredict.exceptions import NotReachable
from passpredict import *

from .utils import assert_datetime_approx

try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo


def assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end, dt_tol=1):
    satellite = SatellitePredictor.from_tle(tle)
    observer = Observer(location, satellite)
    expected_overpasses = observer.pass_list(start, end, method='brute', aos_at_dg=10, time_step=5, tol=0.1)
    observer = Observer(location, satellite)
    date = start
    if len(expected_overpasses) == 0:
        with pytest.raises(NotReachable):
            pass_ = observer.next_pass(date, limit_date=end)
        return
    for expected_pass in expected_overpasses:
        pass_ = observer.next_pass(date, limit_date=end, aos_at_dg=10, tol=0.5)
        assert_datetime_approx(pass_.aos.dt, expected_pass.aos.dt, dt_tol)
        assert_datetime_approx(pass_.los.dt, expected_pass.los.dt, dt_tol)
        assert_datetime_approx(pass_.tca.dt, expected_pass.tca.dt, dt_tol)
        expected_duration = (expected_pass.los.dt - expected_pass.aos.dt).total_seconds()
        assert pass_.duration == pytest.approx(expected_duration, abs=dt_tol*2)
        assert pass_.tca.elevation == pytest.approx(expected_pass.tca.elevation, abs=0.5)
        date = pass_.los.dt + timedelta(minutes=1)


def test_sat_14129_observer_error():
    """
    satid: 14129 TLE: ['1 14129U 83058B   22033.24895832 -.00000141  00000+0  00000+0 0  9997', '2 14129  26.3866 118.1312 5979870  45.7131 349.8958  2.05868779262630']
    Date: 2/2/2022 - 2/3/2022
    Date: 2022-02-08 21:39:13.586360+00:00
    Location: Austin, TX  30.2711, -97.7437
    Error: could not find ascending phase

    Satellite 14129 has a highly eccentric orbit
    """
    tle = TLE(14129, (
        '1 14129U 83058B   22033.24895832 -.00000141  00000+0  00000+0 0  9997',
        '2 14129  26.3866 118.1312 5979870  45.7131 349.8958  2.05868779262630'
    ))
    satellite = SGP4Predictor.from_tle(tle)
    location = Location("Austin, TX", 30.2711, -97.7437, 0)
    observer = Observer(location, satellite)
    start = datetime(2022, 2, 9, 14, 0, tzinfo=ZoneInfo("America/Chicago"))
    with pytest.raises(Exception):
        pass_ = observer.next_pass(start, aos_at_dg=10, tol=0.75)

    # start = datetime(2022, 2, 8, 21, 15, tzinfo=ZoneInfo("America/Chicago"))
    # end = start + timedelta(days=2)
    # observer.pass_list(start, limit_date=end)
