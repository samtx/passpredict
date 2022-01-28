import datetime
from math import floor

import click
from rich.console import Console
from rich.table import Table
from rich.align import Align

from . import __version__
from .sources import CelestrakTLESource
from .geocoding import NominatimGeocoder
from .satellites import SGP4Predictor
from .locations import Location
from .observers import PassType, Observer


@click.command()
@click.option('-s', '--satellite-id', type=int)  # satellite id
@click.option('-d', '--days', default=10, type=int) # day range
@click.option('-loc', '--location', 'location_query', default='', type=str)  # For geocoding location query
@click.option('-lat', '--latitude', default='', type=str)   # latitude
@click.option('-lon', '--longitude', default='', type=str)  # longitude
@click.option('-h', '--height', default=0.0, type=float)    # height
@click.option('--twelve/--twentyfour', '-12/-24', is_flag=True, default=False)  # 12 hour / 24 hour format
@click.option('-a', '--all', 'alltypes', is_flag=True, default=False)  # show all pass types
@click.option('-q', '--quiet', is_flag=True, default=False)
@click.option('-v', '--verbose', is_flag=True, default=False)
@click.option('--summary', is_flag=True, default=False)  # make summary table of results
@click.option('--no-cache', is_flag=True, default=False)
@click.version_option(version=__version__)
def main(satellite_id, days, location_query, latitude, longitude, height, twelve, alltypes, quiet, verbose, summary, no_cache):
    """
    Command line interface for pass predictions
    """
    if location_query:
        location = NominatimGeocoder.query(location_query)
    else:
        location = Location(
            latitude_deg=float(latitude),
            longitude_deg=float(longitude),
            elevation_m=float(height),
            name=""
        )
    satid = satellite_id
    date_start = datetime.datetime.now(tz=location.timezone)
    date_end = date_start + datetime.timedelta(days=days)
    min_elevation = 10.0 # degrees

    source = CelestrakTLESource()
    source.load()
    tle = source.get_tle(satid)
    satellite = SGP4Predictor.from_tle(tle)
    observer = Observer(location, satellite, aos_at_dg=min_elevation, tolerance_s=0.75)
    visible_only = not alltypes
    pass_iterator = observer.iter_passes(date_start, date_end, visible_only=visible_only)
    overpasses = list(pass_iterator)
    # Filter visible passes only unless all passes are requested
    manager = PasspredictManager(
        location,
        satellite,
        tle,
        twelve,
        alltypes,
        quiet,
        verbose,
        summary,
    )
    console = Console()
    if not quiet:
        header = manager.make_results_header()
        console.print(header, highlight=False)
    if len(overpasses) == 0:
        console.print(f"No {'visible' if visible_only else ''} overpasses found")
    else:
        table = manager.overpass_table(overpasses)
        console.print(table)
    source.save()
    return 0


class PasspredictManager:
    """
    Manager to hold all options regarding CLI passpredict query
    """
    def __init__(self,
        location,
        satellite,
        tle,
        twelvehour = False,
        alltypes = False,
        quiet = False,
        verbose = False,
        summary = False,
    ):
        self.location = location
        self.satellite = satellite
        self.tle = tle
        self.twelvehour = twelvehour
        self.alltypes = alltypes
        self.quiet = quiet
        self.verbose = verbose
        self.summary = summary

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
        satellite_id = self.tle.satid
        # Print datetimes with the correct timezone
        line = f"Satellite ID {satellite_id} "
        if self.tle.name:
            line += f"{self.tle.name} "
        line += "overpasses "
        if self.location.name:
            line += f"for {self.location.name:s}"
        line += "\n"
        header += line
        header += f"Lat={self.location.lat:.4f}\u00B0, Lon={self.location.lon:.4f}\u00B0, Timezone {self.location.timezone}\n"
        header += f"Using TLE with epoch {self.tle.epoch.isoformat()}\n"
        header += f"{self.tle.tle1:s}\n"
        header += f"{self.tle.tle2:s}\n"
        return header

    def make_summary_table(self, overpasses):
        """
        Make a summary data table to print to console
        """
        table = Table()
        table.add_column(Align("Date", "center"), justify="right")
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
            row += self.point_string(overpass.aos)
            row += self.point_string(overpass.tca)
            row += self.point_string(overpass.los)
            if overpass.type:
                row.append(overpass.type.value)# + '\n'
                fg = 'green' if overpass.type.value == PassType.visible else None
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


if __name__ == "__main__":
    main()