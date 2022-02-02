import datetime

from passpredict import CelestrakTLESource, Location, SGP4Predictor, Observer

try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo


def celestrak_source():
    location = Location('Austin, TX', 30.2711, -97.7437, 0)
    date_start = datetime.datetime.now(tz=ZoneInfo('America/Chicago'))
    date_end = date_start + datetime.timedelta(days=10)
    source = CelestrakTLESource()
    tle = source.get_tle(25544)  # International space station, Norad ID 25544
    satellite = SGP4Predictor.from_tle(tle)
    observer = Observer(location, satellite)
    overpasses =  observer.pass_list(date_start, limit_date=date_end)
    print(overpasses)


if __name__ == "__main__":
    celestrak_source()