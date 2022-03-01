import datetime

from rich.console import Console

from passpredict import Location, MemoryTLESource, SatellitePredictor, TLE
from passpredict.observers import Observer
from passpredict.cli import PasspredictManager


def brute_force_observer():
    """
    Run brute force predictor for ISS over Austin, Texas
    """
    location = Location('Austin, TX', 30.2711, -97.7437, 0)
    date_start = datetime.datetime(2021, 10, 2, tzinfo=datetime.timezone.utc)
    date_end = date_start + datetime.timedelta(days=10)
    min_elevation = 10.0  # degrees
    tle = TLE(
        25544,  # International space station, Norad ID 25544
        (
            '1 25544U 98067A   21274.89460679  .00005555  00000-0  10931-3 0  9992',
            '2 25544  51.6449 175.1112 0004190  48.8354  53.9444 15.48895782305133'
        )
    )
    source = MemoryTLESource()
    source.add_tle(tle)
    satellite = SatellitePredictor(25544, source)
    observer = Observer(location, satellite)
    overpasses = observer.pass_list(
        date_start, limit_date=date_end, method='brute',
        aos_at_dg=min_elevation, time_step=5, tol=0.1
    )
    manager = PasspredictManager(
        location,
        tle,
        twelvehour=False,
        quiet=False,
        summary=False,
        alltypes=False,
        verbose=False,
    )
    console = Console()
    console.print(manager.make_results_header())
    console.print(manager.overpass_table(overpasses))
    return


if __name__ == "__main__":
    brute_force_observer()
