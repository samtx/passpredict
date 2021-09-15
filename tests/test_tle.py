from datetime import datetime, timezone

from pytest import approx

import passpredict.tle
from passpredict import OMM

def test_Tle():
    """
    From Celestrak

    ISS (ZARYA)
    1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999
    2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675
    """
    satid = 25544
    epoch = datetime(2020, 6, 14, 23, 36, 58, 495104, tzinfo=timezone.utc)
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    tle_model = passpredict.tle.Tle(tle1=tle1, tle2=tle2)
    tle_schema = tle_model.to_schema()
    assert tle_schema.tle1 == tle1
    assert tle_schema.tle2 == tle2
    assert tle_schema.epoch == epoch
    assert tle_schema.satid == satid


def test_Tle_epoch_proptery():
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    epoch = datetime(2020, 6, 14, 23, 36, 58, 495104, tzinfo=timezone.utc)
    tle = passpredict.tle.Tle(tle1, tle2)
    assert tle.epoch == epoch


def test_Tle_satid_proptery():
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    satid = 25544
    tle = passpredict.tle.Tle(tle1, tle2)
    assert tle.satid == satid


def test_OMM_from_tle():
    tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
    tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
    omm = OMM.from_tle(tle1, tle2)
    assert omm.epoch == approx(2459015.48401036)
    assert omm.jdsatepoch == approx(2459014.5)
    assert omm.inclo == approx(51.6466)
    assert omm.ecco == approx(.0002481)
    assert omm.mo == approx(97.0831)
    assert omm.argpo == approx(58.1246)
    assert omm.bstar == approx(1.7092)
    assert omm.elnum == 9999
    assert omm.revnum == 31675
    assert omm.no_kozai == approx(5.494441482)
    assert omm.nodeo == approx(59.3724)
