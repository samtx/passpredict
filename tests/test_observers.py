from datetime import datetime, timedelta, timezone

import pytest

from passpredict import MemoryTLESource, Observer, Location
from passpredict.exceptions import NotReachable
from passpredict import *



def assert_datetime_approx(dt1, dt2, delta_seconds):
    """
    Compare two python datetimes. Assert difference is <= delta_seconds
    """
    diff = (dt1 - dt2).total_seconds()
    assert diff == pytest.approx(0.0, abs=delta_seconds)


def assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end, dt_tol=1):
    satellite = SatellitePredictor.from_tle(tle)
    brute_force_observer = BruteForceObserver(location, satellite, aos_at_dg=10, time_step=5, tolerance_s=0.1)
    expected_overpasses = list(brute_force_observer.iter_passes(start, end))
    observer = Observer(location, satellite, aos_at_dg=10, tolerance_s=0.5)
    date = start
    if len(expected_overpasses) == 0:
        with pytest.raises(NotReachable):
            pass_ = observer.next_pass(date, limit_date=end)
        return
    for expected_pass in expected_overpasses:
        pass_ = observer.next_pass(date, limit_date=end)
        assert_datetime_approx(pass_.aos.dt, expected_pass.aos.dt, dt_tol)
        assert_datetime_approx(pass_.los.dt, expected_pass.los.dt, dt_tol)
        assert_datetime_approx(pass_.tca.dt, expected_pass.tca.dt, dt_tol)
        expected_duration = (expected_pass.los.dt - expected_pass.aos.dt).total_seconds()
        assert pass_.duration == pytest.approx(expected_duration, abs=dt_tol*2)
        assert pass_.tca.elevation == pytest.approx(expected_pass.tca.elevation, abs=0.5)
        date = pass_.los.dt + timedelta(minutes=1)


@pytest.mark.parametrize('observer_class, observer_kwargs', [
    pytest.param(Observer, {'tolerance_s': 0.5}, id='Standard Observer'),
    pytest.param(BruteForceObserver, {'time_step': 5, 'tolerance_s': 0.5}, id='BruteForceObserver'),
])
def test_bugsat_predictions(observer_class, observer_kwargs):
    """
    From orbit-predictor test suite,  newsat 2
    https://github.com/satellogic/orbit-predictor/blob/master/tests/test_accurate_predictor.py
    """
    tle_lines = (
        "1 40014U 14033E   14294.41438078  .00003468  00000-0  34565-3 0  3930",
        "2 40014  97.9781 190.6418 0032692 299.0467  60.7524 14.91878099 18425"
    )
    satellite = SatellitePredictor.from_tle(TLE('Bugsat', tle_lines))

    # ARG in orbit_predictor.locations
    location = Location("ARG", latitude_deg=-31.2884, longitude_deg=-64.2032868, elevation_m=492.96)
    observer = observer_class(location, satellite, **observer_kwargs)

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
                "2014-10-22 20:18:11.921921", '%Y-%m-%d %H:%M:%S.%f').replace(tzinfo=UTC)

        pass_ = observer.next_pass(date)
        assert_datetime_approx(pass_.aos.dt, aos, 1)
        assert_datetime_approx(pass_.los.dt, los, 1)
        assert_datetime_approx(pass_.tca.dt, max_elevation_date, 1)
        assert pass_.duration == pytest.approx(duration_s, abs=2)
        assert pass_.tca.elevation == pytest.approx(max_elev_deg, abs=0.05)


@pytest.mark.slow
@pytest.mark.parametrize('location', [
    pytest.param(Location('Svalbard Satellite Station, Norway', 78.2297, 15.3975, 458), id='Svalbard'),
    pytest.param(Location('Inuvik, Canada', 68.3607, -133.7230, 15), id='Inuvik'),
    pytest.param(Location('Kiev, Ukraine', 50.4501, 30.5234, 179), id='Kiev'),
    pytest.param(Location('Austin, Texas, USA', 32.1234, -97.1234, 0), id='Austin'),
    pytest.param(Location('Quito, Ecuador', 0.1807, -78.4678, 2850), id='Quito'),
    pytest.param(Location('Johannesburg, South Africa', -26.2041, 28.0473, 1753), id='Johannesburg'),
    pytest.param(Location('New Norcia Deep Space Station, Australia', -31.0483, 116.1914, 252), id='New Norica'),
    pytest.param(Location('McMurdo Station, Antarctica', -77.8419, 166.6863, 10), id='McMurdo'),
])
class TestStandardObserver:

    @pytest.mark.parametrize('tle', [
        pytest.param(TLE('ISS', ('1 25544U 98067A   21274.89460679  .00005555  00000-0  10931-3 0  9992','2 25544  51.6449 175.1112 0004190  48.8354  53.9444 15.48895782305133')), id='ISS'),
        pytest.param(TLE('HST', ('1 20580U 90037B   21274.54126351  .00001079  00000-0  54194-4 0  9991','2 20580  28.4698  83.0715 0002784 122.8787 269.0686 15.09757146527044')), id='HST'),
    ])
    def test_leo_common_satellites(self, location, tle):
        """
        Test common LEO satellite overpasses using brute force observer as truth model
        """
        start = datetime(2021, 10, 2, tzinfo=timezone.utc)
        end = datetime(2021, 10, 4, tzinfo=timezone.utc)
        assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end)


    @pytest.mark.parametrize('tle', [
        pytest.param(TLE('Terra', ('1 25994U 99068A   21274.75027490  .00000107  00000-0  33588-4 0  9990','2 25994  98.1678 347.1838 0001354  91.6968  61.5062 14.57132833158940')), id='Terra'),
        pytest.param(TLE('Envisat', ('1 27386U 02009A   21274.71722274  .00000080  00000-0  39878-4 0  9992','2 27386  98.1579 263.3896 0001195  87.0592  86.1444 14.38044645 26342')), id='Envisat'),
        pytest.param(TLE('Aprizesat 5', ('1 37792U 11044E   21274.60878477  .00000190  00000-0  35454-4 0  9996','2 37792  98.3609 136.8622 0059209 108.1425 252.6243 14.75523908544976')), id='Aprizesat 5'),
    ])
    def test_leo_sso_sun_synchronous_satellites(self, location, tle):
        start = datetime(2021, 10, 2, tzinfo=timezone.utc)
        end = datetime(2021, 10, 4, tzinfo=timezone.utc)
        assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end)


    @pytest.mark.parametrize('tle', [
        pytest.param(TLE('OFEQ 11', ('1 41759U 16056A   20133.81295366 0.00000000  00000-0  00000-0 0    06','2 41759 141.8552 279.8922 0006334 285.1363  74.8629 15.34726396    09')), id='OFEQ 11'),
    ])
    def test_leo_retrograde_satellites(self, location, tle):
        start = datetime(2021, 10, 2, tzinfo=timezone.utc)
        end = datetime(2021, 10, 4, tzinfo=timezone.utc)
        assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end)


@pytest.mark.skip(reason="Need to determine how to count geosynchronous overpasses")
@pytest.mark.parametrize('location',[
    pytest.param(Location('Inuvik, Canada', 68.3607, -133.7230, 15), id='Inuvik'),
    pytest.param(Location('Austin, Texas, USA', 32.1234, -97.1234, 0), id='Austin'),
    pytest.param(Location('Quito, Ecuador', 0.1807, -78.4678, 2850), id='Quito'),
])
@pytest.mark.parametrize('tle, start, end', [
    pytest.param(
        TLE('Intelsat 5', ('1 24916U 97046A   21274.50960067  .00000091  00000-0  00000-0 0  9997','2 24916   6.6785  62.7523 0003074 141.8639 212.2992  1.00271606 88373')),
        datetime(2021, 10, 2, tzinfo=timezone.utc),
        datetime(2021, 10, 4, tzinfo=timezone.utc),
        id='Intelsat 5'
    ),
    pytest.param(
        TLE('Galaxy 19', ('1 33376U 08045A   22021.62511755 -.00000127  00000-0  00000-0 0  9999', '2 33376   0.0154  72.0379 0003076 242.5276 294.4103  1.00270478 48785')),
        datetime(2022, 1, 21, tzinfo=timezone.utc),
        datetime(2022, 1, 23, tzinfo=timezone.utc),
        id='Galaxy 19'
    ),
    pytest.param(
        TLE('MUOS 5', ('1 41622U 16041A   22021.62018662 -.00000097  00000-0  00000-0 0  9992', '2 41622   5.9528 302.0771 0196044 233.2720  63.8350  1.00272532 20963')),
        datetime(2022, 1, 21, tzinfo=timezone.utc),
        datetime(2022, 1, 23, tzinfo=timezone.utc),
        id='MUOS 5'
    ),
])
def test_west_geo_geosynchronous_satellites(location, tle, start, end):
    dt_tol = 2
    assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end, dt_tol=dt_tol)


@pytest.mark.skip(reason="Need to determine how to count geosynchronous overpasses")
@pytest.mark.parametrize('location',[
    pytest.param(Location('Perth, Australia', -31.9523, 115.8613, 0), id='Perth'),
    pytest.param(Location('Singapore', 1.3521, 103.8198, 0), id='Singapore'),
    pytest.param(Location('Seoul, South Korea', 37.5665, 126.9780, 38), id='Seoul'),
])
@pytest.mark.parametrize('tle, start, end', [
    pytest.param(
        TLE('BEIDOU 3 IGSO-2', ('1 44337U 19035A   22021.41416855 -.00000102  00000-0  00000-0 0  9994', '2 44337  55.1112 173.0353 0022302 186.9138  25.6402  1.00298568  9574')),
        datetime(2022, 1, 21, tzinfo=timezone.utc),
        datetime(2022, 1, 23, tzinfo=timezone.utc),
        id='BEIDOU 3'
    ),
])
def test_east_geo_geosynchronous_satellites(location, tle, start, end):
    tol = 2  # datetime tolerance in seconds
    assert_overpass_accuracy_with_brute_force_observer(location, tle, start, end, dt_tol=tol)




if __name__ == "__main__":
    pytest.run(__file__)