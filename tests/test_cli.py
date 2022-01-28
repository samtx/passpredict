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
def satellite(tle):
    """
    List of generated overpasses for cli table output testing
    """
    return SGP4Predictor.from_tle(tle)


@pytest.fixture(scope="module")
def overpasses(location, satellite):
    """
    List of generated overpasses for cli table output testing
    """
    observer = Observer(location, satellite)
    start = datetime.datetime(2022, 1, 27)
    end = start + datetime.timedelta(days=4)
    passes = list(observer.iter_passes(start, end))
    return passes


def test_cli():
    runner = CliRunner()
    result = runner.invoke(cli.main,'-s 25544 --location="austin, texas"')
    assert result.exit_code == 0


@pytest.mark.parametrize('twelve', (
    pytest.param(True, id='12hour'),
    pytest.param(False, id='24hour'),
))
class TestCliTable:
    def test_cli_summary_table(self, location, satellite, tle, overpasses, twelve):
        manager = cli.PasspredictManager(location, satellite, tle,
            twelvehour=twelve, summary=True
        )
        manager.make_summary_table(overpasses)

    def test_cli_detail_table(self, location, satellite, tle, overpasses, twelve):
        manager = cli.PasspredictManager(location, satellite, tle,
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