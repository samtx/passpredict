# from datetime import datetime, timezone

# from astrodynamics import TLE

# def test_Tle_from_string():
#     """
#     From Celestrak

#     ISS (ZARYA)
#     1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999
#     2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675
#     """
#     satid = 25544
#     epoch = datetime(2020, 6, 14, 23, 36, 58, 495104, tzinfo=timezone.utc)
#     tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
#     tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
#     tle = Tle.from_string(tle1=tle1, tle2=tle2)
#     assert tle.tle1 == tle1
#     assert tle.tle2 == tle2
#     assert tle.epoch == epoch
#     assert tle.satid == satid


# def test_Tle_epoch_proptery():
#     tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
#     tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
#     epoch = datetime(2020, 6, 14, 23, 36, 58, 495104, tzinfo=timezone.utc)
#     tle = Tle.from_string(tle1, tle2)
#     assert tle.epoch == epoch


# def test_Tle_satid_proptery():
#     tle1 = '1 25544U 98067A   20166.98401036  .00000505  00000-0  17092-4 0  9999'
#     tle2 = '2 25544  51.6466 359.3724 0002481  58.1246  97.0831 15.49444148231675'
#     satid = 25544
#     tle = Tle.from_string(tle1, tle2)
#     assert tle.satid == satid