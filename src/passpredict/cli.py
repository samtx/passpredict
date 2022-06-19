from __future__ import annotations
import datetime
from math import floor
from typing import Union, Sequence
import warnings

import click
from rich.console import Console
from rich.table import Table
from rich.align import Align

from . import __version__
from .caches import JsonCache
from .exceptions import CelestrakError, PassAlgorithmError, PropagationError
from .sources import CelestrakTLESource
from .geocoding import NominatimGeocoder
from .satellites import SGP4Propagator
from .locations import Location
from .observers import PassType, Observer, Visibility
from .orbit import TLE


@click.command(no_args_is_help=True)
@click.option('-s', '--satid', 'satids', type=int, default=list(), multiple=True, help="Satellite NORAD ID")  # satellite id
@click.option('-c', '--category', 'categories', type=str, default=list(), multiple=True, help="Celestrak category string. Example: 'visual', 'amateur'")  # Celestrak category
@click.option('-d', '--days', default=1, type=int, help="Number of days to compute overpasses") # day range
@click.option('-loc', '--location', 'location_query', default='', type=str, help="Location name string")  # For geocoding location query
@click.option('-lat', '--latitude', default='', type=str, help="Location latitude in decimal degrees North")   # latitude
@click.option('-lon', '--longitude', default='', type=str, help="Location longitude in decimal degrees East")  # longitude
@click.option('-h', '--height', default=0.0, type=float, help="Location height in meters")    # height
@click.option('--min-elevation', default=10, type=click.FloatRange(min=-90, max=90))  # minimum elevation in degrees
@click.option('--twelve/--twentyfour', '-12/-24', is_flag=True, default=False)  # 12 hour / 24 hour format
@click.option('-a', '--all', 'alltypes', is_flag=True, default=False)  # show all pass types
@click.option('-q', '--quiet', is_flag=True, default=False)
@click.option('-v', '--verbose', is_flag=True, default=False)
@click.option('--summary', is_flag=True, default=False)  # make summary table of results
@click.option('--no-cache', is_flag=True, default=False)
@click.version_option(version=__version__)
def main(satids, categories, days, location_query, latitude, longitude, height, min_elevation, twelve, alltypes, quiet, verbose, summary, no_cache):
    """
    Predict satellite passes over a location on Earth.
    """
    cache = JsonCache()
    cache.load()

    if location_query:
        geocoder = NominatimGeocoder(cache=cache)
        location = geocoder.query(location_query)
    elif latitude != '' and longitude != '':
        location = Location(
            latitude_deg=float(latitude),
            longitude_deg=float(longitude),
            elevation_m=float(height),
            name=""
        )
    else:
        raise click.BadParameter("Must specify observing location name with --location or coordinates with --lat and --lon")
    date_start = datetime.datetime.now(tz=location.timezone)
    date_end = date_start + datetime.timedelta(days=days)

    source = CelestrakTLESource(cache=cache)
    # Get TLE data for selected satellites and categories
    tles = []
    if len(satids) == 0 and len(categories) == 0:
        raise click.BadParameter("Must specify satellite with --satid or --category")
    try:
        for satid in satids:
            tles.append(source.get_tle(satid))
        for category in categories:
            tles += source.get_tle_category(category)
    except CelestrakError as err:
        raise click.BadParameter(str(err))

    # Get overpasses for each TLE
    console = Console()

    visible_only = not alltypes
    overpasses = []

    if not quiet:
        msg = f"Computing overpasses for {len(tles)} satellites over {days} days..."
        console.print(msg)

    for tle in tles:
        satellite = SGP4Propagator.from_tle(tle)
        observer = Observer(location, satellite)
        try:
            overpasses += observer.pass_list(
                date_start, date_end, visible_only=visible_only,
                aos_at_dg=min_elevation, tol=0.75
            )
        except PropagationError:
            warnings.warn(f"Propagation Error raised for {tle}. Overpasses skipped.")
        except PassAlgorithmError:
            warnings.warn(f"Pass Algorithm Error raised for {tle}. Overpasses skipped.")


    # Sort overpass list by AOS date
    overpasses.sort(key=lambda pass_: pass_.aos.dt)

    manager = PasspredictManager(
        location,
        tles,
        twelve,
        alltypes,
        quiet,
        verbose,
        summary,
    )
    if not quiet:
        header = manager.make_results_header()
        console.print(header, highlight=False)
    if len(overpasses) == 0:
        console.print(f"No {'visible' if visible_only else ''} overpasses found")
    else:
        table = manager.overpass_table(overpasses)
        console.print(table)
    cache.save()
    return 0


class PasspredictManager:
    """
    Manager to hold all options regarding CLI passpredict query
    """
    def __init__(
        self,
        location,
        tles: Union[TLE, Sequence[TLE]],
        twelvehour=False,
        alltypes=False,
        quiet=False,
        verbose=False,
        summary=False,
    ):
        self.location = location
        self.tles = tles if is_tle_sequence(tles) else [tles]
        self.twelvehour = twelvehour
        self.alltypes = alltypes
        self.quiet = quiet
        self.verbose = verbose
        self.summary = summary
        self.sat_name_idx = {t.satid: t.name for t in self.tles}

    def point_string(self, point):
        time = point.dt.astimezone(self.location.tz)
        point_data = []
        if self.twelvehour:
            point_data.append(time.strftime("%I:%M:%S").lstrip("0") + ' ' + time.strftime("%p").lower())
        else:
            point_data.append(time.strftime("%H:%M:%S"))
        point_data.append("{:>2}\u00B0".format(int(point.elevation)))
        point_data.append("{:3}".format(point.direction))
        return point_data

    def overpass_table(self, overpasses):
        """
        Return a formatted string for tabular output

        Params:
            overpasses: list
                A list of Overpass objects
        Return:
            table : str
                tabular formatted string
        """
        # Data results
        if self.summary:
            table = self.make_summary_table(overpasses)
        else:
            table = self.make_detail_table(overpasses)
        return table

    def make_results_header(self):
        header = ""
        line = ""
        if len(self.tles) == 1:
            tle = self.tles[0]
            # Print datetimes with the correct timezone
            line = f"Satellite ID {tle.satid} "
            if tle.name:
                line += f"{tle.name} "
            line += "overpasses "
        if self.location.name:
            line += f"for {self.location.name:s}"
        line += "\n"
        header += line
        header += f"Lat={self.location.lat:.4f}\u00B0, Lon={self.location.lon:.4f}\u00B0, Timezone {self.location.timezone}\n"
        if len(self.tles) == 1:
            header += f"Using TLE with epoch {tle.epoch.isoformat()}\n"
            header += f"{tle.tle1:s}\n"
            header += f"{tle.tle2:s}\n"
        return header

    def make_summary_table(self, overpasses):
        """
        Make a summary data table to print to console
        """
        table = Table()
        table.add_column(Align("Date", "center"), justify="right")
        table.add_column(Align("Satellite", 'center'), justify="left")
        table.add_column(Align("Duration", "center"), justify="right")
        table.add_column(Align("Max Elev", "center"), justify="right")
        table.add_column(Align("Type", "center"), justify='center')
        for overpass in overpasses:
            row = []
            date = overpass.aos.dt.astimezone(self.location.timezone)
            day = date.strftime("%x").lstrip('0')
            # round to nearest minute
            if date.second >= 30:
                date += datetime.timedelta(minutes=1)
            if self.twelvehour:
                time = date.strftime("%I:%M").lstrip("0")
                ampm = date.strftime('%p').lower()
                row.append(f"{day:>8}  {time:>5} {ampm}")
            else:
                time = date.strftime("%H:%M")
                row.append(f"{day}  {time}")
            row.append(f"{overpass.satid} {self.sat_name_idx[overpass.satid]}")
            row.append(get_min_sec_string(overpass.duration))
            row.append(f"{int(overpass.tca.elevation):2}\u00B0")
            if overpass.type:
                row.append(overpass.type.value)# + '\n'
                fg = 'green' if overpass.type.value == PassType.visible else None
            else:
                fg = None
            table.add_row(*row, style=fg )
        return table

    def make_detail_table(self, overpasses):
        """
        Make a detailed data table to print to console
        """
        table = Table()
        table.add_column(Align("Date", 'center'), justify="right")
        table.add_column(Align("Satellite", 'center'), justify="left")
        for x in ('Start', 'Max', 'End'):
            table.add_column(Align(f"{x}\nTime", "center"), justify="right")
            table.add_column(Align(f"{x}\nEl", "center"), justify="right", width=4)
            table.add_column(Align(f"{x}\nAz", "center"), justify="right", width=4)
        table.add_column(Align("Type", "center"), justify='center')
        for overpass in overpasses:
            row = []
            row.append(overpass.aos.dt.astimezone(self.location.tz).strftime("%x").lstrip('0'))
            # if overpass.brightness is not None:
            #     brightness_str = f"{overpass.brightness:4.1f}"
            # else:
            #     brightness_str = " "*4
            # table_data += " "*2 + brightness_str
            row.append(f"{overpass.satid} {self.sat_name_idx[overpass.satid]}")
            row += self.point_string(overpass.aos)
            row += self.point_string(overpass.tca)
            row += self.point_string(overpass.los)
            if overpass.type:
                row.append(overpass.type.value)# + '\n'
                fg = 'green' if overpass.type.value == Visibility.visible else None
            else:
                fg = None
            table.add_row(*row, style=fg )
        return table


@click.command()
def flush():
    """
    Remove all expired keys from cache
    """
    raise NotImplementedError


def get_min_sec_string(total_seconds: int) -> str:
    """
    Get total number of seconds, return string with min:sec format
    """
    nmin = floor(total_seconds / 60)
    nsec = total_seconds - nmin * 60
    return f"{nmin:.0f}:{nsec:02.0f}"


def is_tle_sequence(x) -> bool:
    """  Return true if x is a list, tuple, or set of TLE objects """
    if isinstance(x, TLE):
        return False
    if isinstance(x, list) and isinstance(x[0], TLE):
        return True
    if isinstance(x, tuple) and isinstance(x[0], TLE):
        return True
    if isinstance(x, set):
        y = tuple(x)
        if isinstance(y[0], TLE):
            return True
        return False
    return False



if __name__ == "__main__":
    main()