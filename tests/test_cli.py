import datetime

import pytest
from click.testing import CliRunner

from passpredict import cli
from passpredict import Location, SGP4Predictor, Observer, TLE



@pytest.fixture(scope="module")
def location():
    return Location("Austin, Texas", 30.2711, -97.7437, 0)


@pytest.fixture(scope="module")
def tle():
    """
    List of generated overpasses for cli table output testing
    """
    tle = TLE(25544, (
        "1 25544U 98067A   22027.75547926  .00007563  00000-0  14194-3 0  9991",
        "2 25544  51.6451 312.0370 0006903  67.9238  24.2396 15.49669141323407"
    ), name="ISS")
    return tle


@pytest.fixture(scope="module")
def overpasses(location, tle):
    """
    List of generated overpasses for cli table output testing
    """
    satellite = SGP4Predictor.from_tle(tle)
    observer = Observer(location, satellite)
    start = datetime.datetime(2022, 1, 27)
    end = start + datetime.timedelta(days=4)
    passes = list(observer.iter_passes(start, end))
    return passes

@pytest.mark.parametrize('options_string',[
    pytest.param('-s 25544 --location="austin, texas"', id="one sat, geocode location"),
    pytest.param('-s 25544 -lat 30.1234 -lon -98.134', id="one sat, location coordinates"),
    pytest.param('--category visual --location="austin, texas" -d 1', id='sat category'),
    pytest.param('-s 25544 -s 20580 -s 27386 --location="austin, texas" -d 3', id="multi sat"),
])
def test_cli(options_string):
    runner = CliRunner()
    result = runner.invoke(cli.main, options_string)
    assert result.exit_code == 0


@pytest.mark.parametrize('twelve', (
    pytest.param(True, id='12hour'),
    pytest.param(False, id='24hour'),
))
class TestCliTable:
    def test_cli_summary_table(self, location, tle, overpasses, twelve):
        manager = cli.PasspredictManager(location, tle,
            twelvehour=twelve, summary=True
        )
        manager.make_summary_table(overpasses)

    def test_cli_detail_table(self, location, tle, overpasses, twelve):
        manager = cli.PasspredictManager(location, tle,
            twelvehour=twelve, summary=False
        )
        manager.make_detail_table(overpasses)


@pytest.mark.parametrize('nsec, expected_str', [
    (60, '1:00'),
    (59, '0:59'),
    (121, '2:01'),
    (306, '5:06'),
])
def test_get_min_sec_string(nsec, expected_str):
    s = cli.get_min_sec_string(nsec)
    assert s == expected_str


@pytest.mark.parametrize('tles, res', [
    pytest.param([TLE(123, ('123', '123'))], True, id='list of one'),
    pytest.param([TLE(123, ('123', '123')), TLE(456, ('456', '456'))], True, id='list of two'),
    pytest.param((TLE(123, ('123', '123')), TLE(456, ('456', '456'))), True, id='tuple of two'),
    pytest.param((TLE(123, ('123', '123')),), True, id='tuple of one'),
    pytest.param({TLE(123, ('123', '123')), TLE(456, ('456', '456'))}, True, id='set of two'),
    pytest.param(TLE(123, ('123', '123')), False, id='single TLE'),
    pytest.param('blahblah', False, id='string, not TLE'),
])
def test_is_tle_sequence(tles, res):
    assert cli.is_tle_sequence(tles) == res