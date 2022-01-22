from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import pytest

from passpredict import Observer, Location, SGP4Predictor
from passpredict import *



def assert_datetime_approx(dt1, dt2, delta_seconds):
    """
    Compare two python datetimes. Assert difference is <= delta_seconds
    """
    diff = (dt1 - dt2).total_seconds()
    assert diff == pytest.approx(0.0, abs=delta_seconds)


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
        # tca = datetime.strptime(f"{date_str} {line_parts[5]}+0100", "%Y %d %b %H:%M:%S%z")
        vis_end = datetime.strptime(f"{date_str} {line_parts[8]}+0100", "%Y %d %b %H:%M:%S%z")

        # Find next pass
        pass_ = observer.get_next_pass(date, limit_date=end, visible=True)
        assert_datetime_approx(pass_.vis_begin.dt, vis_begin, 1.5)
        assert_datetime_approx(pass_.vis_end.dt, vis_end, 1.5)
        # assert_datetime_approx(pass_.tca.dt, tca, 5)
        date = pass_.los.dt + timedelta(minutes=10)


if __name__ == "__main__":
    pytest.main([__file__])