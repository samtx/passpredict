import datetime
from zoneinfo import ZoneInfo

from passpredict import CelestrakTLESource, Location, SGP4Predictor, Observer


def celestrak_source():
    location = Location('Austin, TX', 30.2711, -97.7437, 0)
    date_start = datetime.datetime.now(tz=ZoneInfo('America/Chicago'))
    date_end = date_start + datetime.timedelta(days=10)
    source = CelestrakTLESource()
    tle = source.get_tle(25544)  # International space station, Norad ID 25544
    satellite = SGP4Predictor.from_tle(tle)
    observer = Observer(location, satellite)
    pass_iterator =  observer.iter_passes(date_start, limit_date=date_end)
    overpasses = list(pass_iterator)
    print(overpasses)


if __name__ == "__main__":
    celestrak_source()