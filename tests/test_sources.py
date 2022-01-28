
import pytest

from passpredict import sources
from passpredict import TLE


class TestCelestrakSource:

    def test_celestrak_source_get_tle(self, tmp_path):
        cache_path = tmp_path / 'tle.json'
        source = sources.CelestrakTLESource(cache_path)
        tle = source.get_tle(25544)
        assert tle.satid == 25544
        assert len(tle.lines) == 2
        assert tle.name == 'ISS (ZARYA)'
        assert isinstance(tle, TLE)
