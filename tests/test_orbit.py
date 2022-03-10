import pytest
from pytest import approx
from sgp4.api import Satrec, WGS84

from passpredict import SGP4Propagator, Orbit, TLE
from passpredict.sources import parse_multiple_tles


tles = """DELTA 2 R/B(1)
1 20453U 90008B   22062.60993960  .00004513  00000+0  24453-3 0  9998
2 20453  35.6228 170.5595 0271539 122.7895 239.9285 14.87880488704123
COSMOS 2058
1 20465U 90010A   22062.38664240  .00001287  00000+0  10517-3 0  9996
2 20465  82.4896  32.7736 0013613 256.3753 103.5955 14.98536050721749
SL-14 R/B
1 20466U 90010B   22062.45103841  .00000221  00000+0  25019-4 0  9991
2 20466  82.5045 202.7650 0021943  87.8471 272.5255 14.81484304711116
SL-14 R/B
1 20511U 90018B   22062.51964734  .00000195  00000+0  22678-4 0  9991
2 20511  82.5254 151.8946 0016337 227.6694 142.7939 14.80135863705429
HST
1 20580U 90037B   22062.55741664  .00001715  00000+0  91794-4 0  9991
2 20580  28.4694 150.2239 0002454 341.5292  39.8862 15.10085625550192
SL-16 R/B
1 20625U 90046B   22062.58369218  .00000064  00000+0  58003-4 0  9997
2 20625  70.9997 138.1527 0013998   4.2620 355.8619 14.14687418640886
COSMOS 2084
1 20663U 90055A   22062.59591716  .00000576  00000+0  68150-4 0  9992
2 20663  62.7920 243.5960 0106378  38.7301 322.1357 14.92522736716492
INTELSAT 4-F3
1 05709U 71116A   22063.50151007  .00000005  00000+0  00000+0 0  9992
2 05709  10.1102 306.7060 0013109  82.1726 147.8397  0.99624918103304
SL-6 R/B(2)
1 20666U 90055D   22062.61018495  .00000578  00000+0  65727-4 0  9995
2 20666  62.7842 158.3620 0124294  61.9786 299.3816 14.92792395718459
MOLNIYA 2-9
1 07276U 74026A   22061.50516359  .00000037  00000+0  00000+0 0  9994
2 07276  64.2662 242.3347 6535050 282.6627  15.9855  2.45095606246598
MOLNIYA 2-10
1 07376U 74056A   22061.00843876 -.00000824  00000+0  00000+0 0  9995
2 07376  63.9571 105.8642 6701222 287.6291  13.4463  2.01100127348197
SWISSCUBE
1 35932U 09051B   22062.87666685  .00000386  00000+0  97300-4 0  9990
2 35932  98.5818 275.7932 0007949  19.8655 340.2850 14.56704149660725
BEESAT-1
1 35933U 09051C   22063.52050935  .00000369  00000+0  93424-4 0  9993
2 35933  98.5763 278.1522 0006340  32.1833 327.9750 14.56741763660978"""


def generate_tle_params(tle_strings):
    tles = parse_multiple_tles(tle_strings)
    for tle in tles:
        yield pytest.param(tle, id=tle.name)


@pytest.mark.parametrize('tle', generate_tle_params(tles.splitlines()))
def test_satellite_from_tle(tle):
    sat = SGP4Propagator.from_tle(tle)._propagator
    sat2 = Satrec.twoline2rv(tle.tle1, tle.tle2, WGS84)
    assert sat.bstar == approx(sat2.bstar)
    assert sat.argpo == approx(sat2.argpo)
    assert sat.ndot == approx(sat2.ndot)
    assert sat.nddot == approx(sat2.nddot)
    assert sat.mo == approx(sat2.mo)
    assert sat.nodeo == approx(sat2.nodeo)
    assert sat.no_kozai == approx(sat2.no_kozai)
    assert sat.inclo == approx(sat2.inclo)
    assert sat.jdsatepoch == approx(sat2.jdsatepoch)
    assert sat.jdsatepochF == approx(sat2.jdsatepochF)







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