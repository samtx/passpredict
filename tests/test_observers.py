from datetime import datetime, timedelta, timezone
import json
from collections import defaultdict
import pathlib
import itertools

import pytest

from passpredict import MemoryTLESource, Observer
from passpredict import *
from passpredict.observers.predicted_pass import PassPoint
from .data.generate_predictions import generate_predictions

# From orbit-predictor test suite
# https://github.com/satellogic/orbit-predictor/blob/master/tests/test_accurate_predictor.py
SAT_ID = '41558U'  # newsat 2
LINES = (
    '1 41558U 16033C   17065.21129769  .00002236  00000-0  88307-4 0  9995',
    '2 41558  97.4729 144.7611 0014207  16.2820 343.8872 15.26500433 42718',
)

ONE_SECOND = timedelta(seconds=1)

def assert_datetime_approx(dt1, dt2, delta_seconds):
    """
    Compare two python datetimes. Assert difference is <= delta_seconds
    """
    diff = (dt1 - dt2).total_seconds()
    assert diff == pytest.approx(0.0, abs=delta_seconds)


def test_satid_bugsat_predictions():
    """
    From orbit-predictor test suite,  newsat 2
    https://github.com/satellogic/orbit-predictor/blob/master/tests/test_accurate_predictor.py
    """
    start = datetime(2017, 3, 6, 7, 51)
    satid = ' -1'
    tle_lines = (
        "1 40014U 14033E   14294.41438078  .00003468  00000-0  34565-3 0  3930",
        "2 40014  97.9781 190.6418 0032692 299.0467  60.7524 14.91878099 18425"
    )
    # tle = TLE(satid, tle_lines, start)
    source = MemoryTLESource()
    source.add_tle(TLE(satid, tle_lines))
    satellite = SatellitePredictor(satid, source)

    # ARG in orbit_predictor.locations
    location = Location("ARG", latitude_deg=-31.2884, longitude_deg=-64.2032868, elevation_m=492.96)

    observer = Observer(location, satellite)

    STK_DATA = """
    ------------------------------------------------------------------------------------------------
        AOS                      TCA                      LOS                      Duration      Max El
    ------------------------------------------------------------------------------------------------
        2014/10/23 01:27:33.224  2014/10/23 01:32:41.074  2014/10/23 01:37:47.944  00:10:14.720   12.76
        2014/10/23 03:01:37.007  2014/10/23 03:07:48.890  2014/10/23 03:14:01.451  00:12:24.000   39.32
        2014/10/23 14:49:34.783  2014/10/23 14:55:44.394  2014/10/23 15:01:51.154  00:12:16.000   41.75
        2014/10/23 16:25:54.939  2014/10/23 16:30:50.152  2014/10/23 16:35:44.984  00:09:50.000   11.45
        2014/10/24 01:35:47.889  2014/10/24 01:41:13.181  2014/10/24 01:46:37.548  00:10:50.000   16.07
        2014/10/24 03:10:23.486  2014/10/24 03:16:27.230  2014/10/24 03:22:31.865  00:12:08.000   30.62
        2014/10/24 14:58:07.378  2014/10/24 15:04:21.721  2014/10/24 15:10:33.546  00:12:26.000   54.83
        2014/10/24 16:34:48.635  2014/10/24 16:39:20.960  2014/10/24 16:43:53.204  00:09:04.000    8.78
        2014/10/25 01:44:05.771  2014/10/25 01:49:45.487  2014/10/25 01:55:24.414  00:11:18.000   20.07
        2014/10/25 03:19:12.611  2014/10/25 03:25:05.674  2014/10/25 03:30:59.815  00:11:47.000   24.09"""  # NOQA

    for line in STK_DATA.splitlines()[4:]:
        line_parts = line.split()
        aos = datetime.strptime(" ".join(line_parts[:2]), '%Y/%m/%d %H:%M:%S.%f')
        max_elevation_date = datetime.strptime(" ".join(line_parts[2:4]),
                                                    '%Y/%m/%d %H:%M:%S.%f')
        los = datetime.strptime(" ".join(line_parts[4:6]), '%Y/%m/%d %H:%M:%S.%f')
        for d in (aos, max_elevation_date, los):
            d = d.replace(tzinfo=timezone.utc)
        duration = datetime.strptime(line_parts[6], '%H:%M:%S.%f')
        duration_s = timedelta(
            minutes=duration.minute, seconds=duration.second).total_seconds()
        max_elev_deg = float(line_parts[7])

        try:
            date = pass_.los.dt  # NOQA
        except UnboundLocalError:
            date = datetime.strptime(
                "2014-10-22 20:18:11.921921", '%Y-%m-%d %H:%M:%S.%f')

        pass_ = observer.get_next_pass(date)
        assert_datetime_approx(pass_.aos.dt, aos, 1)
        assert_datetime_approx(pass_.los.dt, los, 1)
        assert_datetime_approx(pass_.tca.dt, max_elevation_date, 1)
        assert pass_.duration == pytest.approx(duration_s, abs=2)
        assert pass_.tca.elevation == pytest.approx(max_elev_deg, abs=0.05)


def test_satid_bugsat_brute_force_predictions():
    """
    From orbit-predictor test suite,  newsat 2
    https://github.com/satellogic/orbit-predictor/blob/master/tests/test_accurate_predictor.py
    """
    start = datetime(2017, 3, 6, 7, 51)
    satid = ' -1'
    tle_lines = (
        "1 40014U 14033E   14294.41438078  .00003468  00000-0  34565-3 0  3930",
        "2 40014  97.9781 190.6418 0032692 299.0467  60.7524 14.91878099 18425"
    )
    # tle = TLE(satid, tle_lines, start)
    source = MemoryTLESource()
    source.add_tle(TLE(satid, tle_lines))
    satellite = SatellitePredictor(satid, source)

    # ARG in orbit_predictor.locations
    location = Location("ARG", latitude_deg=-31.2884, longitude_deg=-64.2032868, elevation_m=492.96)

    observer = BruteForceObserver(location, satellite, time_step=5, tolerance_s=0.1)

    STK_DATA = """
    ------------------------------------------------------------------------------------------------
        AOS                      TCA                      LOS                      Duration      Max El
    ------------------------------------------------------------------------------------------------
        2014/10/23 01:27:33.224  2014/10/23 01:32:41.074  2014/10/23 01:37:47.944  00:10:14.720   12.76
        2014/10/23 03:01:37.007  2014/10/23 03:07:48.890  2014/10/23 03:14:01.451  00:12:24.000   39.32
        2014/10/23 14:49:34.783  2014/10/23 14:55:44.394  2014/10/23 15:01:51.154  00:12:16.000   41.75
        2014/10/23 16:25:54.939  2014/10/23 16:30:50.152  2014/10/23 16:35:44.984  00:09:50.000   11.45
        2014/10/24 01:35:47.889  2014/10/24 01:41:13.181  2014/10/24 01:46:37.548  00:10:50.000   16.07
        2014/10/24 03:10:23.486  2014/10/24 03:16:27.230  2014/10/24 03:22:31.865  00:12:08.000   30.62
        2014/10/24 14:58:07.378  2014/10/24 15:04:21.721  2014/10/24 15:10:33.546  00:12:26.000   54.83
        2014/10/24 16:34:48.635  2014/10/24 16:39:20.960  2014/10/24 16:43:53.204  00:09:04.000    8.78
        2014/10/25 01:44:05.771  2014/10/25 01:49:45.487  2014/10/25 01:55:24.414  00:11:18.000   20.07
        2014/10/25 03:19:12.611  2014/10/25 03:25:05.674  2014/10/25 03:30:59.815  00:11:47.000   24.09"""  # NOQA

    UTC = timezone.utc
    for line in STK_DATA.splitlines()[4:]:
        line_parts = line.split()
        aos = datetime.strptime(" ".join(line_parts[:2]), '%Y/%m/%d %H:%M:%S.%f').replace(tzinfo=UTC)
        max_elevation_date = datetime.strptime(" ".join(line_parts[2:4]),
                                                    '%Y/%m/%d %H:%M:%S.%f').replace(tzinfo=UTC)
        los = datetime.strptime(" ".join(line_parts[4:6]), '%Y/%m/%d %H:%M:%S.%f').replace(tzinfo=UTC)
        duration = datetime.strptime(line_parts[6], '%H:%M:%S.%f')
        duration_s = timedelta(
            minutes=duration.minute, seconds=duration.second).total_seconds()
        max_elev_deg = float(line_parts[7])

        try:
            date = pass_.los.dt + timedelta(minutes=1) # NOQA
        except UnboundLocalError:
            date = datetime.strptime(
                "2014-10-22 20:18:11.921921", '%Y-%m-%d %H:%M:%S.%f')

        pass_ = observer.get_next_pass(date)
        assert_datetime_approx(pass_.aos.dt, aos, 0.5)
        assert_datetime_approx(pass_.los.dt, los, 0.5)
        assert_datetime_approx(pass_.tca.dt, max_elevation_date, 0.5)
        assert pass_.duration == pytest.approx(duration_s, abs=1)
        assert pass_.tca.elevation == pytest.approx(max_elev_deg, abs=0.05)


def generate_prediction_parameters():
    """
    Generate test observations using BruteForceObserver
    """
    # Get predictions.json if it exists
    fname = 'predictions.json'
    fpath = pathlib.Path(__file__).parent / 'data' / fname
    if not fpath.exists():
        generate_predictions(fname)
    else:
        with open(fpath, 'r') as f:
            data = json.load(f)

    # Encode locations list
    locations = []
    for l in data['locations']:
        location = Location(l['name'], l['lat'], l['lon'], l['h'])
        locations.append(location)

    # Encode satellites list
    satellites = []
    for s in data['satellites']:
        source = MemoryTLESource()
        source.add_tle(TLE(s['satid'], s['lines']))
        satellite = SatellitePredictor(s['satid'], source)
        satellite.name = s['name']
        satellites.append(satellite)

    # Encode overpasses dictionary
    overpasses = defaultdict(list)
    for ls_key in data['overpasses']:
        for op in data['overpasses'][ls_key]:
            aos = PassPoint(
                datetime.fromisoformat(op['aos']['dt']),
                float(op['aos']['range']),
                float(op['aos']['azimuth']),
                float(op['aos']['elevation'])
            )
            tca = PassPoint(
                datetime.fromisoformat(op['tca']['dt']),
                float(op['tca']['range']),
                float(op['tca']['azimuth']),
                float(op['tca']['elevation'])
            )
            los = PassPoint(
                datetime.fromisoformat(op['los']['dt']),
                float(op['los']['range']),
                float(op['los']['azimuth']),
                float(op['los']['elevation'])
            )
            pass_ = PredictedPass(
                op['satid'],
                op['location'],
                aos,
                tca,
                los,
            )
            overpasses[ls_key].append(pass_)

    start = datetime.fromisoformat(data['start'])
    params = []
    for l, s in itertools.product(locations, satellites):
        key = f"{l.name}-{s.name}"

        # Currently, the Inuvik, Canada-Intelsat 5 combination is throwing an error.
        # Mark as xfail for now
        if key == "Inuvik, Canada-Intelsat 5":
            param = pytest.param(l, s, overpasses[key], start, marks=pytest.mark.xfail, id=key)
        else:
            param = pytest.param(l, s, overpasses[key], start, id=key)
        params.append(param)
    return params


@pytest.mark.parametrize(
    'location, satellite, overpasses, start',
    generate_prediction_parameters()
)
def test_observer_with_prediction_suite(location, satellite, overpasses, start):
    date = start
    observer = Observer(location, satellite, aos_at_dg=10, tolerance_s=0.5)
    for expected_pass in overpasses:
        pass_ = observer.get_next_pass(date)
        assert_datetime_approx(pass_.aos.dt, expected_pass.aos.dt, 1)
        assert_datetime_approx(pass_.los.dt, expected_pass.los.dt, 1)
        assert_datetime_approx(pass_.tca.dt, expected_pass.tca.dt, 1)
        expected_duration = (expected_pass.los.dt - expected_pass.aos.dt).total_seconds()
        assert pass_.duration == pytest.approx(expected_duration, abs=2)
        assert pass_.tca.elevation == pytest.approx(expected_pass.tca.elevation, abs=0.5)
        date = pass_.los.dt + timedelta(minutes=1)


if __name__ == "__main__":
    pytest.run(__file__)


