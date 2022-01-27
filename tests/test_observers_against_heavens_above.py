from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import pytest
from pytest import approx

from passpredict import Observer, Location, SGP4Predictor
from passpredict import *
from passpredict.observers.base import PassPoint



def assert_datetime_approx(dt1, dt2, delta_seconds):
    """
    Compare two python datetimes. Assert difference is <= delta_seconds
    """
    diff = (dt1 - dt2).total_seconds()
    assert diff == approx(0.0, abs=delta_seconds)


def assert_passpoint_approx(pt1, pt2, *, dt_tol=1, el_tol=1, az_tol=10, range_tol=10):
    """
    Compare two PassPoint objects
    """
    assert_datetime_approx(pt1.dt, pt2.dt, dt_tol)
    assert pt1.elevation == approx(pt2.elevation, abs=el_tol)
    assert pt1.azimuth == approx(pt2.azimuth, abs=az_tol)
    assert pt1.direction == pt2.direction
    assert pt1.range == approx(pt2.range, abs=range_tol)


def test_heavens_above_zurich_iss_visibility_predictions():
    """
    From heavens-above.com. Queried 1/21/2022
    """
    tle_lines = (
        "1 25544U 98067A   22021.68747065 -.00000457  00000-0  00000+0 0  9996",
        "2 25544  51.6270 342.0904 0006847  44.2808  14.3181 15.49598315322461"
    )
    satellite = SGP4Predictor.from_tle(TLE('ISS', tle_lines))

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
        pass_ = observer.get_next_pass(date, limit_date=end, visible_only=True)
        assert_datetime_approx(pass_.vis_begin.dt, vis_begin, 1.5)
        assert_datetime_approx(pass_.vis_end.dt, vis_end, 1.5)
        assert_datetime_approx(pass_.vis_tca.dt, tca, 1.5)
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
            vis_begin = datetime.strptime(f"{date_str} {line_parts[2]}-0600", "%Y %d %b %H:%M:%S%z")
            tca = datetime.strptime(f"{date_str} {line_parts[5]}-0600", "%Y %d %b %H:%M:%S%z")
            vis_end = datetime.strptime(f"{date_str} {line_parts[8]}-0600", "%Y %d %b %H:%M:%S%z")

            # Find next pass
            pass_ = self.observer.get_next_pass(date, limit_date=end, visible_only=True)
            assert_datetime_approx(pass_.vis_begin.dt, vis_begin, 1.5)
            assert_datetime_approx(pass_.vis_end.dt, vis_end, 1.5)
            assert_datetime_approx(pass_.vis_tca.dt, tca, 1.5)
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
        # aos_pt.sun_elevation = -5.9
        tca_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 42, 3, tzinfo=TZ),
            range=891,
            azimuth=161,
            elevation=34,
            brightness=1.9
        )
        # tca_pt.sun_elevation = -7.2
        los_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 45, 46, tzinfo=TZ),
            range=1775,
            azimuth=96,
            elevation=10,
            brightness=2.8
        )
        # los_pt.sun_elevation = -8.0
        vis_end_pt = PassPoint(
            dt=datetime(2022, 1, 31, 18, 46, 3, tzinfo=TZ),
            range=1876,
            azimuth=95,
            elevation=9,
            brightness=2.9
        )
        # vis_end_pt.sun_elevation = -8.0

        # Find next pass
        pass_ = self.observer.get_next_pass(start, visible_only=True)
        tol = {
            'dt_tol': 1.5,
            'el_tol': 1,
            'az_tol': 10,
            'range_tol': 10,
        }
        assert_passpoint_approx(pass_.aos, aos_pt, **tol)
        assert_passpoint_approx(pass_.tca, tca_pt, **tol)
        assert_passpoint_approx(pass_.los, los_pt, **tol)



if __name__ == "__main__":
    pytest.main([__file__])