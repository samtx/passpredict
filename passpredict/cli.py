import datetime
from math import floor

import click
from rich.console import Console
from rich.table import Table
from rich.align import Align

from .caches import JsonCache
from .tle import get_TLE, Tle
from .core import predict_single_satellite_overpasses
from ._time import julian_date
from .satellites import SatellitePredictor
from .locations import Location
from .observers import PassType
from .sources import TLE


@click.command()
@click.option('-s', '--satellite-id', type=int)  # satellite id
@click.option('-d', '--days', default=10, type=click.IntRange(1, 14, clamp=True)) # day range
@click.option('-lat', '--latitude', type=float)  # latitude
@click.option('-lon', '--longitude', type=float)  # longitude
@click.option('-h', '--height', default=0.0, type=float)    # height
@click.option('--twelve/--twentyfour', '-12/-24', is_flag=True, default=False)  # 12 hour / 24 hour format
@click.option('-a', '--all', 'alltypes', is_flag=True, default=False)  # show all pass types
@click.option('-q', '--quiet', is_flag=True, default=False)
@click.option('-v', '--verbose', is_flag=True, default=False)
@click.option('--summary', is_flag=True, default=False)  # make summary table of results
@click.option('--no-cache', is_flag=True, default=False)
def main(satellite_id, days, latitude, longitude, height, twelve, alltypes, quiet, verbose, summary, no_cache):
    """
    Command line interface for pass predictions
    """
    location = Location(
        latitude_deg=latitude,
        longitude_deg=longitude,
        elevation_m=height,
        name=""
    )
    satid = satellite_id
    date_start = datetime.datetime.now(tz=datetime.timezone.utc)
    min_elevation = 10.0 # degrees

    if no_cache:
        tle = get_TLE(satid)
    else:
        with JsonCache() as cache:
            sat_key = f'sat:{satid}'
            res = cache.get(sat_key)
            if res:
                tle = Tle(tle1=res['tle1'], tle2=res['tle2'])
            else:
                tle = get_TLE(satid)
                cache.set(sat_key, tle.dict(), ttl=86400)


    satellite = SatellitePredictor(satid)
    satellite.tle = TLE(tle.satid, (tle.tle1, tle.tle2), tle.epoch)
    satellite.set_propagator()
    overpasses = predict_single_satellite_overpasses(satellite, location, date_start, days, min_elevation=min_elevation)
    # Filter visible passes only unless all passes are requested
    if not alltypes:
        overpasses = filter(lambda p: p.type == PassType.visible, overpasses)
    overpass_table(overpasses, location, tle, tz=location.timezone, twelvehour=twelve, quiet=quiet, summary=summary)
    return 0


def overpass_table(overpasses, location, tle, tz=None, twelvehour=False, quiet=False, summary=False):
    """
    Return a formatted string for tabular output

    Params:
        overpasses: list
            A list of Overpass objects
        location: Location object
        tle: Tle object

    Return:
        table : str
            tabular formatted string
"""
    table_title = ""
    if not quiet:
        satellite_id = tle.satid
        # Print datetimes with the correct timezone
        table_title += f"Satellite ID {satellite_id} overpasses for {location.name:s}\n"
        table_title += f"Lat={location.lat:.4f}\u00B0, Lon={location.lon:.4f}\u00B0, Timezone {tz}\n"
        table_title += f"Using TLE with epoch {tle.epoch.isoformat()}\n"
        table_title += f"{tle.tle1:s}\n"
        table_title += f"{tle.tle2:s}\n"
    click.secho(table_title)

    # Data results
    if summary:
        table = make_summary_table(overpasses, tz, twelvehour)
    else:
        table = make_detail_table(overpasses, tz, twelvehour)
    console = Console()
    console.print(table)


def make_summary_table(overpasses, tz, twelvehour):
    """
    Make a summary data table to print to console
    """
    table = Table()
    table.add_column(Align("Date", "center"), justify="right")
    table.add_column(Align("Duration", "center"), justify="right")
    table.add_column(Align("Max Elev", "center"), justify="right")
    table.add_column(Align("Type", "center"), justify='center')

    def get_min_sec_string(total_seconds: int) -> str:
        """
        Get total number of seconds, return string with min:sec format
        """
        nmin = floor(total_seconds / 60)
        nsec = total_seconds - nmin * 60
        return f"{nmin:.0f}:{nsec:0.0f}"

    for overpass in overpasses:
        row = []
        date = overpass.aos.dt.astimezone(tz)
        day = date.strftime("%x").lstrip('0')
        # round to nearest minute
        if date.second >= 30:
            date.replace(minute=date.minute+1)
        if twelvehour:
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


def make_detail_table(overpasses, tz, twelvehour):
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

    def point_string(point):
        time = point.dt.astimezone(tz)
        point_data = []
        if twelvehour:
            point_data.append(time.strftime("%I:%M:%S").lstrip("0") + ' ' + time.strftime("%p").lower())
        else:
            point_data.append(time.strftime("%H:%M:%S"))
        point_data.append("{:>2}\u00B0".format(int(point.elevation)))
        point_data.append("{:3}".format(point.direction))
        return point_data

    for overpass in overpasses:
        row = []
        row.append(overpass.aos.dt.astimezone(tz).strftime("%x").lstrip('0'))
        # if overpass.brightness is not None:
        #     brightness_str = f"{overpass.brightness:4.1f}"
        # else:
        #     brightness_str = " "*4
        # table_data += " "*2 + brightness_str
        row += point_string(overpass.aos)
        row += point_string(overpass.tca)
        row += point_string(overpass.los)
        if overpass.type:
            row.append(overpass.type.value)# + '\n'
            fg = 'green' if overpass.type.value == PassType.visible else None
        else:
            fg = None
        table.add_row(*row, style=fg )
    return table


def echo(*a, **kw):
    if 'end' in kw:
        kw.pop('end')
        kw.update({'nl': False})
    return click.echo(*a, **kw)


@click.command()
def flush():
    """
    Remove all expired keys from cache
    """
    pass

if __name__ == "__main__":
    main()