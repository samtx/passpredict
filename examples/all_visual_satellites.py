import datetime
try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo

from passpredict import CelestrakTLESource, Location, SGP4Predictor, Observer


def all_visual_satellites():
    location = Location('Austin, TX', 30.2711, -97.7437, 0)
    source = CelestrakTLESource()
    tles = source.get_tle_category("visual")  # Visible satellites
    overpasses = []
    start = datetime.datetime.now(tz=ZoneInfo('America/Chicago'))
    end = start + datetime.timedelta(days=1)
    for tle in tles:
        satellite = SGP4Predictor.from_tle(tle)
        observer = Observer(location, satellite)
        passes = observer.pass_list(start, end, visible_only=True, aos_at_dg=10)
        overpasses += passes
    return overpasses


if __name__ == "__main__":
    from rich.pretty import pprint
    overpasses = all_visual_satellites()
    pprint(overpasses)