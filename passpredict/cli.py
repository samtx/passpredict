import datetime
from typing import NamedTuple

import click

from .schemas import Satellite, Location, PassType
from .geocoding import geocoder
from .utils import Cache
from .tle import get_TLE
from .predictions import predict

@click.command()
@click.option('-s', '--satellite-id', type=int)  # satellite id
@click.option('-l', '--location', 'location_string', default="", type=str)  # location string
@click.option('-u', '--utc-offset', default=0.0, type=float, prompt=True)  # utc offset
@click.option('-d', '--days', default=10, type=click.IntRange(1, 14, clamp=True)) # day range
@click.option('-lat', '--latitude', type=float)  # latitude
@click.option('-lon', '--longitude', type=float)  # longitude
@click.option('-h', '--height', default=0.0, type=float)    # height
@click.option('--twelve/--twentyfour', '-12/-24', is_flag=True, default=False)  # 12 hour / 24 hour format
@click.option('-a', '--all', 'alltypes', is_flag=True, default=False)  # show all pass types
@click.option('-q', '--quiet', is_flag=True, default=False)
@click.option('-v', '--verbose', is_flag=True, default=False)
@click.option('--no-cache', is_flag=True, default=False)
def main(satellite_id, location_string, utc_offset, days, latitude, longitude, height, twelve, alltypes, quiet, verbose, no_cache):
    """
    Command line interface for pass predictions
    """
    if latitude and longitude:
        lat = round(float(latitude), 4)
        lon = round(float(longitude), 4)
    else:
        data = geocoder(location_string)      # Prompt for location
        lat = round(float(data['lat']), 4)
        lon = round(float(data['lon']), 4)

    tz = datetime.timezone(datetime.timedelta(hours=utc_offset))
    location = Location(
        lat=lat,
        lon=lon,
        h=height,
        name=location_string
    )
    satellite = Satellite(
        id=satellite_id,
    )
    tle = get_TLE(satellite.id)
    date_start = datetime.date.today()
    date_end = date_start + datetime.timedelta(days=days)
    min_elevation = 10.01 # degrees
    cache = Cache() if not no_cache else None
    visible_only = not alltypes
    overpasses = predict(location, satellite, date_start=date_start, date_end=date_end, dt_seconds=1, min_elevation=min_elevation, verbose=verbose, cache=cache, print_fn=echo, visible_only=visible_only)
    overpass_table(overpasses, location, tle, tz, twelvehour=twelve, quiet=quiet)
    if cache is not None:
        with cache:
            cache.flush()
            
    return 0


def overpass_table(overpasses, location, tle, tz=None, twelvehour=False, quiet=False):
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
    table_header = ""
    if not quiet:
        satellite_id = tle.satid
        # Print datetimes with the correct timezone
        table_title += f"Satellite ID {satellite_id} overpasses for {location.name:s}\n"
        table_title += f"Lat={location.lat:.4f}\u00B0, Lon={location.lon:.4f}\u00B0, Timezone {tz}\n"
        table_title += f"Using TLE with epoch {tle.epoch.isoformat()}\n"
        table_title += f"{tle.tle1:s}\n"
        table_title += f"{tle.tle2:s}\n"
    click.secho(table_title)
    if twelvehour:
        point_header = "  Time    El  Az "
        point_header_underline = "--------- --- ---"
    else:
        #   Time   Elx Azx
        point_header = "  Time   El  Az "
        point_header_underline = "-------- --- ---"
    table_header +=  f"           {'Start':^17s}   {'Maximum':^17s}   {'End':^17s}\n"
    table_header +=  "  Date    Mag  {0}   {0}   {0}      Type\n".format(point_header)
    table_header += "--------  ----  "
    table_header += point_header_underline + " "*3
    table_header += point_header_underline + " "*3
    table_header += point_header_underline + " "*3
    table_header += "-"*10
    click.echo(table_header)

    def point_string(point):
        time = point.datetime.astimezone(tz)
        if twelvehour:
            point_line = time.strftime("%I:%M:%S") + time.strftime("%p")[0].lower()
        else:
            point_line = time.strftime("%H:%M:%S")
        point_line += " " + "{:>2}\u00B0".format(int(point.elevation))
        point_line += " " + "{:3}".format(point.direction_from_azimuth())
        return point_line

    for overpass in overpasses:
        table_data = "{}".format(overpass.start_pt.datetime.astimezone(tz).strftime("%m/%d/%y"))
        if overpass.brightness is not None:
            brightness_str = f"{overpass.brightness:4.1f}"
        else:
            brightness_str = " "*4
        table_data += " "*2 + brightness_str
        table_data += " "*2 + point_string(overpass.start_pt) + ' |'
        table_data += " " + point_string(overpass.max_pt) + ' |'
        table_data += " " + point_string(overpass.end_pt)
        table_data += " "*4 + f'{overpass.type.value:^9}'# + '\n'
        fg = 'green' if overpass.type.value == PassType.visible else None
        click.secho(table_data, fg=fg)

    
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