from datetime import datetime, timedelta, timezone

import pytest
from pytest import approx

from passpredict import Observer, Location, SGP4Predictor
from passpredict import *
from passpredict.observers import PassPoint, Visibility
from passpredict.zoneinfo import ZoneInfo



def assert_datetime_approx(dt1, dt2, delta_seconds):
    """
    Compare two python datetimes. Assert difference is <= delta_seconds
    """
    diff = (dt1 - dt2).total_seconds()
    assert diff == approx(0.0, abs=delta_seconds)


def assert_passpoint_approx(pt1, pt2, *, dt_tol=1, el_tol=1, az_tol=10, range_tol=10, bright_tol=3):
    """
    Compare two PassPoint objects
    """
    assert_datetime_approx(pt1.dt, pt2.dt, dt_tol)
    assert pt1.elevation == approx(pt2.elevation, abs=el_tol)
    assert pt1.azimuth == approx(pt2.azimuth, abs=az_tol)
    assert pt1.direction == pt2.direction
    assert pt1.range == approx(pt2.range, abs=range_tol)
    # if (pt1.brightness is not None) or (pt2.brightness is not None):
    #     assert pt1.brightness == approx(pt2.brightness, abs=bright_tol)


def assert_sun_elevation_at_date(location, dt, el, tol=0.1):
    """
    Compare the sun elevation at the location
    """
    assert location.sun_elevation(dt.astimezone(timezone.utc)) == approx(el, abs=tol)

def test_heavens_above_zurich_iss_visibility_predictions():
    """
    From heavens-above.com. Queried 1/21/2022
    """
    tle_lines = (
        "1 25544U 98067A   22021.68747065 -.00000457  00000-0  00000+0 0  9996",
        "2 25544  51.6270 342.0904 0006847  44.2808  14.3181 15.49598315322461"
    )
    satellite = SGP4Predictor.from_tle(TLE('ISS', tle_lines))
    satellite.intrinsic_mag = -1.8

    location = Location("Zurich", 47.3744, 8.5410, 0)
    observer = Observer(location, satellite, aos_at_dg=10, tolerance_s=0.5)

    HEAVENS_ABOVE_PREDICTIONS = """
    Date	Brightness	Start	Highest point	End	Pass type
    (mag)	Time	Alt.	Az.	Time	Alt.	Az.	Time	Alt.	Az.
    21 Jan	-4.0	19:02:20	10°	WSW	19:05:41	88°	NNW	19:05:54	76°	ENE	visible
    22 Jan	-3.6	18:14:24	10°	SW	18:17:41	64°	SSE	18:19:57	19°	ENE	visible
    22 Jan	-1.4	19:51:23	10°	W	19:52:53	23°	WNW	19:52:53	23°	WNW	visible
    23 Jan	-3.4	19:03:14	10°	W	19:06:31	51°	NNW	19:06:53	48°	NNE	visible
    24 Jan	-3.6	18:15:08	10°	WSW	18:18:27	66°	NNW	18:20:52	17°	ENE	visible
    24 Jan	-1.3	19:52:24	10°	WNW	19:53:48	21°	WNW	19:53:48	21°	WNW	visible
    25 Jan	-3.0	19:04:17	10°	WNW	19:07:27	39°	N	19:07:46	38°	N	visible
    26 Jan	-3.1	18:16:07	10°	W	18:19:21	43°	N	18:21:43	16°	ENE	visible
    26 Jan	-1.3	19:53:19	10°	WNW	19:54:39	21°	NW	19:54:39	21°	NW	visible
    27 Jan	-3.0	19:05:15	10°	WNW	19:08:25	38°	N	19:08:38	38°	NNE	visible
    28 Jan	-2.9	18:17:10	10°	WNW	18:20:19	37°	N	18:22:38	16°	ENE	visible
    28 Jan	-1.5	19:54:06	10°	WNW	19:55:34	23°	WNW	19:55:34	23°	WNW	visible
    29 Jan	-3.5	19:06:05	10°	WNW	19:09:22	50°	NNE	19:09:36	48°	NNE	visible
    30 Jan	-3.2	18:18:04	10°	WNW	18:21:17	42°	N	18:23:42	16°	ENE	visible
    30 Jan	-1.7	19:54:52	10°	WNW	19:56:38	28°	WNW	19:56:38	28°	WNW	visible"""  # NOQA

    TZ = ZoneInfo('Europe/Zurich')
    UTC = timezone.utc
    start = datetime(2022, 1, 21, tzinfo=TZ).astimezone(UTC)
    end = datetime(2022, 1, 31, tzinfo=TZ).astimezone(UTC)
    date = start
    for line in HEAVENS_ABOVE_PREDICTIONS.splitlines()[3:]:
        line_parts = line.split('\t')
        date_str = f"2022 {line_parts[0]}"
        vis_begin = datetime.strptime(f"{date_str} {line_parts[2]}+0100", "%Y %d %b %H:%M:%S%z")
        tca = datetime.strptime(f"{date_str} {line_parts[5]}+0100", "%Y %d %b %H:%M:%S%z")
        vis_end = datetime.strptime(f"{date_str} {line_parts[8]}+0100", "%Y %d %b %H:%M:%S%z")

        # Find next pass
        pass_ = observer.next_pass(date, limit_date=end, visible_only=True)
        assert_datetime_approx(pass_.vis_begin.dt, vis_begin, 1.5)
        assert_datetime_approx(pass_.vis_end.dt, vis_end, 1.5)
        assert_datetime_approx(pass_.vis_tca.dt, tca, 1.5)
        assert isinstance(pass_.brightness, float)
        assert pass_.type == Visibility.visible
        date = pass_.los.dt + timedelta(minutes=10)


@pytest.fixture(scope='function')
def sanantonio_location():
    """
    Location object for San Antonio, Texas
    """
    return Location("San Antonio, Texas", 29.4246, -98.4951, 0)


@pytest.fixture(scope='function')
def hst_tle():
    """
    TLE for Hubble Space Telescope
    From heavens-above.com. Queried 1/24/2022
    """
    tle_lines = (
        "1 20580U 90037B   22024.50840557  .00001221  00000-0  62432-4 0  9995",
        "2 20580  28.4712  42.1300 0002330 285.0338 208.0692 15.09965616544703"
    )
    return TLE(20580, tle_lines, name="HST")


class TestHeavensAboveSanAntonioHSTVisibilty:

    @classmethod
    def setup_class(cls):
        cls.location = Location("San Antonio, Texas", 29.4246, -98.4951, 0)
        tle_lines = (
            "1 20580U 90037B   22024.50840557  .00001221  00000-0  62432-4 0  9995",
            "2 20580  28.4712  42.1300 0002330 285.0338 208.0692 15.09965616544703"
        )
        tle = TLE(20580, tle_lines, name="HST")
        cls.satellite = SGP4Predictor.from_tle(tle)
        cls.satellite.intrinsic_mag = 2.2
        cls.observer = Observer(cls.location, cls.satellite, aos_at_dg=10, tolerance_s=0.5)

    def test_visibile_overpasses(self):
        """
        From heavens-above.com. Queried 1/24/2022
        Location: San Antonio, Texas, USA
        Satellite: Hubble Space Telescope, ID 20580
        """
        HEAVENS_ABOVE_PREDICTIONS = """
        Date	Brightness	Start	Highest point	End	Pass type
        (mag)	Time	Alt.	Az.	Time	Alt.	Az.	Time	Alt.	Az.
        26 Jan	3.2	19:37:23	10°	S	19:38:32	11°	SSE	19:38:32	11°	SSE	visible
        27 Jan	2.8	19:25:03	10°	S	19:27:23	15°	SSE	19:28:05	14°	SE	visible
        28 Jan	2.6	19:13:09	10°	SSW	19:16:03	18°	SSE	19:17:37	15°	SE	visible
        28 Jan	4.0	20:52:50	10°	WSW	20:52:55	10°	WSW	20:52:55	10°	WSW	visible
        29 Jan	2.4	19:01:26	10°	SSW	19:04:42	23°	SSE	19:07:07	14°	ESE	visible
        29 Jan	3.5	20:41:30	10°	WSW	20:42:25	16°	WSW	20:42:25	16°	WSW	visible
        30 Jan	2.1	18:49:50	10°	SW	18:53:22	28°	SSE	18:56:36	12°	ESE	visible
        30 Jan	3.0	20:30:12	10°	WSW	20:31:54	24°	WSW	20:31:54	24°	WSW	visible
        31 Jan	1.9	18:38:19	10°	SW	18:42:03	34°	SSE	18:45:47	10°	E	visible
        31 Jan	2.3	20:18:54	10°	WSW	20:21:21	34°	WSW	20:21:21	34°	WSW	visible
        01 Feb	1.5	20:07:38	10°	W	20:10:47	51°	WSW	20:10:47	51°	WSW	visible
        02 Feb	0.8	19:56:22	10°	W	20:00:13	74°	SW	20:00:13	74°	SW	visible"""  # NOQA

        TZ = self.location.timezone
        UTC = timezone.utc
        start = datetime(2022, 1, 24, tzinfo=TZ).astimezone(UTC)
        end = datetime(2022, 2, 3, tzinfo=TZ).astimezone(UTC)
        date = start
        for line in HEAVENS_ABOVE_PREDICTIONS.splitlines()[3:]:
            line_parts = line.split('\t')
            date_str = f"2022 {line_parts[0]}"
            # brightness = float(line_parts[1])
            vis_begin = datetime.strptime(f"{date_str} {line_parts[2]}-0600", "%Y %d %b %H:%M:%S%z")
            tca = datetime.strptime(f"{date_str} {line_parts[5]}-0600", "%Y %d %b %H:%M:%S%z")
            vis_end = datetime.strptime(f"{date_str} {line_parts[8]}-0600", "%Y %d %b %H:%M:%S%z")

            # Find next pass
            pass_ = self.observer.next_pass(date, limit_date=end, visible_only=True)
            assert_datetime_approx(pass_.vis_begin.dt, vis_begin, 1.5)
            assert_datetime_approx(pass_.vis_end.dt, vis_end, 1.5)
            assert_datetime_approx(pass_.vis_tca.dt, tca, 1.5)
            assert isinstance(pass_.brightness, float)
            assert pass_.type == Visibility.visible
            # assert pass_.brightness == approx(brightness, abs=3)  # this is very inaccurate
            date = pass_.los.dt + timedelta(minutes=10)


    def test_visibility_one_overpass_Jan_31_2022_1835(self):
        HEAVENS_ABOVE_PREDICTIONS = """
        Date:	31 January 2022
        Orbit:	535 x 538 km, 28.5° (Epoch: 24 January)
        Event	Time	Altitude	Azimuth	Distance (km)	Brightness	Sun altitude
        Rises	18:35:57	0°	236° (SW)	2,658	6.0	-5.9°
        Reaches altitude 10°	18:38:20	10°	226° (SW)	1,770	4.6	-6.4°
        Maximum altitude	18:42:03	34°	161° (SSE)	891	1.9	-7.2°
        Drops below altitude 10°	18:45:46	10°	96° (E)	1,775	2.8	-8.0°
        Enters shadow	18:46:03	9°	95° (E)	1,876	2.9	-8.0°"""
        TZ = self.location.timezone
        UTC = timezone.utc
        start = datetime(2022, 1, 31, 18, 0, 0, tzinfo=TZ).astimezone(UTC)

        aos_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 38, 20, tzinfo=TZ),
            range=1770,
            azimuth=236,
            elevation=10,
            brightness=6.0
        )
        assert_sun_elevation_at_date(self.location, aos_pt.dt, -6.4)  # aos_pt.sun_elevation = -6.4
        tca_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 42, 3, tzinfo=TZ),
            range=891,
            azimuth=161,
            elevation=34,
            brightness=1.9
        )
        assert_sun_elevation_at_date(self.location, tca_pt.dt, -7.2)  # tca_pt.sun_elevation = -7.2
        los_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 45, 46, tzinfo=TZ),
            range=1775,
            azimuth=96,
            elevation=10,
            brightness=2.8
        )
        assert_sun_elevation_at_date(self.location, los_pt.dt, -8.0)  # los_pt.sun_elevation = -8.0
        vis_end_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 46, 3, tzinfo=TZ),
            range=1876,
            azimuth=95,
            elevation=9,
            brightness=2.9
        )
        assert_sun_elevation_at_date(self.location, vis_end_pt.dt, -8.0)   # vis_end_pt.sun_elevation = -8.0

        # Find next pass
        pass_ = self.observer.next_pass(start, visible_only=True)
        tol = {
            'dt_tol': 1.5,
            'el_tol': 1,
            'az_tol': 10,
            'range_tol': 10,
            'bright_tol': 3,
        }
        assert_passpoint_approx(pass_.aos, aos_pt, **tol)
        assert_passpoint_approx(pass_.tca, tca_pt, **tol)
        assert_passpoint_approx(pass_.los, los_pt, **tol)
        assert_passpoint_approx(pass_.vis_end, los_pt, **tol)
        assert pass_.type == Visibility.visible


class TestHeavensAboveCapeTownEnvisatVisibilty:

    @classmethod
    def setup_class(cls):
        cls.location = Location("Cape Town, South Africa", -33.9290, 18.4174, 0)
        tle_lines = (
            "1 27386U 02009A   22030.56291447  .00000091  00000-0  43706-4 0  9996",
            "2 27386  98.1695  18.2467 0001188  87.4527  20.8530 14.38072512 43711"
        )
        tle = TLE(27386, tle_lines, name="Envisat")
        cls.satellite = SGP4Predictor.from_tle(tle)
        cls.satellite.intrinsic_mag = 3.7
        cls.observer = Observer(cls.location, cls.satellite, aos_at_dg=10, tolerance_s=0.5, sunrise_dg=-2)

    def test_visibile_overpasses(self):
        """
        From heavens-above.com. Queried 1/30/2022
        Location: Cape Town, South Africa
        Satellite: Envisat, ID 27386

        Search period start:	31 January 2022 00:00
        Search period end:	8 February 2022 00:00
        Orbit:	764 x 766 km, 98.2° (Epoch: 30 January)

        """
        HEAVENS_ABOVE_PREDICTIONS = """
        Date	Brightness	Start	Highest point	End	Pass type
        (mag)	Time	Alt.	Az.	Time	Alt.	Az.	Time	Alt.	Az.
        31 Jan	5.8	03:39:46	10°	ENE	03:42:47	16°	ESE	03:45:49	10°	SSE	visible
        31 Jan	3.2	05:16:21	10°	N	05:21:17	52°	W	05:26:19	10°	SSW	visible
        01 Feb	3.4	04:39:24	10°	NNE	04:44:26	66°	ESE	04:49:33	10°	S	visible
        02 Feb	5.0	04:03:30	10°	NE	04:07:47	28°	ESE	04:12:07	10°	SSE	visible
        02 Feb	3.9	05:42:19	10°	NNW	05:46:40	28°	W	05:51:03	10°	SW	visible
        03 Feb	6.0	03:29:20	10°	E	03:31:19	12°	ESE	03:33:19	10°	SE	visible
        03 Feb	3.0	05:04:37	10°	N	05:09:41	68°	W	05:14:49	10°	SSW	visible
        04 Feb	3.9	04:27:58	10°	NNE	04:32:53	50°	ESE	04:37:51	10°	S	visible
        05 Feb	5.4	03:52:27	10°	ENE	03:56:17	22°	ESE	04:00:10	10°	SSE	visible
        05 Feb	3.6	05:30:20	10°	NNW	05:35:01	36°	W	05:39:45	10°	SSW	visible
        06 Feb	3.0	04:52:59	10°	NNE	04:58:05	88°	WNW	05:03:15	10°	SSW	visible
        07 Feb	4.4	04:16:39	10°	NE	04:21:21	39°	ESE	04:26:06	10°	S	visible"""  # NOQA

        TZ = self.location.timezone
        UTC = timezone.utc
        start = datetime(2022, 1, 31, 3, tzinfo=TZ)
        end = datetime(2022, 2, 8, tzinfo=TZ)
        date = start
        for line in HEAVENS_ABOVE_PREDICTIONS.splitlines()[3:]:
            line_parts = line.split('\t')
            date_str = f"2022 {line_parts[0]}"
            # brightness = float(line_parts[1])
            vis_begin = datetime.strptime(f"{date_str} {line_parts[2]}+0200", "%Y %d %b %H:%M:%S%z")
            tca = datetime.strptime(f"{date_str} {line_parts[5]}+0200", "%Y %d %b %H:%M:%S%z")
            vis_end = datetime.strptime(f"{date_str} {line_parts[8]}+0200", "%Y %d %b %H:%M:%S%z")

            # Find next pass
            pass_ = self.observer.next_pass(date, limit_date=end, visible_only=True)
            assert_datetime_approx(pass_.vis_begin.dt, vis_begin, 1.5)
            assert_datetime_approx(pass_.vis_end.dt, vis_end, 1.5)
            assert_datetime_approx(pass_.vis_tca.dt, tca, 1.5)
            assert isinstance(pass_.brightness, float)
            assert pass_.type == Visibility.visible
            # assert pass_.brightness == approx(brightness, abs=3)  # this is very inaccurate
            date = pass_.los.dt + timedelta(minutes=10)


    def test_visibility_one_overpass_Feb_06_2022_0450(self):
        HEAVENS_ABOVE_PREDICTIONS = """
        Date:	06 February 2022
        Orbit:	764 x 766 km, 98.2° (Epoch: 30 January)
        Event	Time	Altitude	Azimuth	Distance (km)	Brightness	Sun altitude
        Rises	04:50:42	0°	13° (NNE)	3,219	6.0	-16.1°
        Reaches altitude 10°	04:52:59	10°	13° (NNE)	2,301	5.3	-15.7°
        Maximum altitude	04:58:05	88°	284° (WNW)	779	3.0	-14.8°
        Drops below altitude 10°	05:03:15	10°	193° (SSW)	2,332	5.7	-13.9°
        Sets	05:05:34	0°	192° (SSW)	3,266	6.5	-13.5°"""
        TZ = self.location.timezone
        UTC = timezone.utc
        start = datetime(2022, 2, 6, 0, 0, 0, tzinfo=TZ).astimezone(UTC)

        aos_pt = PassPoint(
            dt=datetime(2022, 2, 6, 4, 52, 59, tzinfo=TZ),
            range=2301,
            azimuth=13,
            elevation=10,
            brightness=5.3
        )
        assert_sun_elevation_at_date(self.location, aos_pt.dt, -15.7, 0.1)  # sun_elevation = -15.7
        tca_pt = PassPoint(
            dt=datetime(2022, 2, 6, 4, 58, 5, tzinfo=TZ),
            range=779,
            azimuth=284,
            elevation=88,
            brightness=3.0
        )
        assert_sun_elevation_at_date(self.location, tca_pt.dt, -14.8)   # tca_pt.sun_elevation = -14.8
        los_pt = PassPoint(
            dt=datetime(2022, 2, 6, 5, 3, 15, tzinfo=TZ),
            range=2332,
            azimuth=193,
            elevation=10,
            brightness=5.7
        )
        assert_sun_elevation_at_date(self.location, los_pt.dt, -13.9)  # los_pt.sun_elevation = -13.9

        # Find next pass
        pass_ = self.observer.next_pass(start, visible_only=True)
        tol = {
            'dt_tol': 1.5,
            'el_tol': 1,
            'az_tol': 10,
            'range_tol': 10,
            'bright_tol': 3,
        }
        assert_passpoint_approx(pass_.vis_begin, aos_pt, **tol)
        assert_passpoint_approx(pass_.vis_tca, tca_pt, **tol)
        assert_passpoint_approx(pass_.vis_end, los_pt, **tol)
        assert pass_.type == Visibility.visible

    def test_visibility_one_overpass_Feb_02_2022_0539(self):
        """
        Date:	02 February 2022
        Orbit:	764 x 766 km, 98.2° (Epoch: 30 January)
        Event	Time	Altitude	Azimuth	Distance (km)	Brightness	Sun altitude
        Rises	05:39:45	0°	347° (NNW)	3,219	5.9	-6.4°
        Reaches altitude 10°	05:42:20	10°	335° (NNW)	2,302	5.0	-5.9°
        Maximum altitude	05:46:40	28°	276° (W)	1,419	3.9	-5.1°
        Drops below altitude 10°	05:51:03	10°	218° (SW)	2,329	5.4	-4.2°
        Sets	05:53:41	0°	206° (SSW)	3,264	6.2	-3.7°"""
        TZ = self.location.timezone
        UTC = timezone.utc
        start = datetime(2022, 2, 2, 4, 30, 0, tzinfo=TZ)  # just after a previous pass

        aos_pt = PassPoint(
            dt=datetime(2022, 2, 2, 5, 42, 20, tzinfo=TZ),
            range=2302,
            azimuth=335,
            elevation=10,
            brightness=5.0
        )
        assert_sun_elevation_at_date(self.location, aos_pt.dt, -5.9)  # aos_pt.sun_elevation = -5.9
        tca_pt = PassPoint(
            dt=datetime(2022, 2, 2, 5, 46, 40, tzinfo=TZ),
            range=1419,
            azimuth=276,
            elevation=28,
            brightness=3.9
        )
        assert_sun_elevation_at_date(self.location, tca_pt.dt, -5.1)  # tca_pt.sun_elevation = -5.1
        los_pt = PassPoint(
            dt=datetime(2022, 2, 2, 5, 51, 3, tzinfo=TZ),
            range=2329,
            azimuth=218,
            elevation=10,
            brightness=5.4
        )
        assert_sun_elevation_at_date(self.location, los_pt.dt, -4.2)  # los_pt.sun_elevation = -4.2

        # Find next pass
        pass_ = self.observer.next_pass(start, visible_only=True)
        tol = {
            'dt_tol': 1.5,
            'el_tol': 1,
            'az_tol': 10,
            'range_tol': 10,
            'bright_tol': 3,
        }
        assert_passpoint_approx(pass_.vis_begin, aos_pt, **tol)
        assert_passpoint_approx(pass_.vis_tca, tca_pt, **tol)
        assert_passpoint_approx(pass_.vis_end, los_pt, **tol)
        assert pass_.type == Visibility.visible


if __name__ == "__main__":
    pytest.main([__file__])