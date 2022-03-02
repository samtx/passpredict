# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import datetime
import cProfile
import pstats

from passpredict import Location, Observer, SGP4Propagator, TLE


# Setup objects
start = datetime.datetime(2020, 7, 14, tzinfo=datetime.timezone.utc)
end = start + datetime.timedelta(days=10)
satid = 'ISS'
tle_lines = (
    "1 25544U 98067A   20196.51422950 -.00000046  00000-0  72206-5 0  9999",
    "2 25544  51.6443 213.2207 0001423 114.8006 342.8278 15.49514729236251"
)
tle = TLE(satid, tle_lines)
min_elevation = 10
satellite = SGP4Propagator.from_tle(tle)
location = Location("Austin, Texas", 30.2711, -97.7434, 0)
observer = Observer(location, satellite)


class PredictOverpasses:
    """
    ISS prediction over Austin, Texas
    """

    def time_observer_iter_passes(self):
        observer.pass_list(start, limit_date=end, aos_at_dg=min_elevation)

    def track_elevation_at_function_calls(self):
        with cProfile.Profile() as pr:
            observer.pass_list(start, limit_date=end, aos_at_dg=min_elevation)
        stats = pstats.Stats(pr)
        fname = '_elevation_mjd'
        keys = stats.stats.keys()
        for k in keys:
            if k[2] == fname:
                break
        ncalls = stats.stats[k][0]
        return ncalls

    def track_elevation_at_function_cache_ratio(self):
        observer.pass_list(start, limit_date=end, aos_at_dg=min_elevation)
        res = observer._elevation_mjd.cache_info()
        return res.hits / (res.misses + res.hits)


class PredictOverpassesBruteForce:
    """
    ISS prediction over Austin, Texas with Brute Force algorithm
    """
    params = (
        [40, 30, 20, 10, 5],
        [1.0, 0.25, 0.1],
    )
    param_names = ['time_step', 'tol']

    def time_brute_force_observer(self, time_step, tol):
        observer.pass_list(
            start, limit_date=end, method='brute',
            aos_at_dg=min_elevation,
            time_step=time_step,
            tol=tol,
        )
