import datetime

import click

from .schemas import Satellite, PassType
from .caches import JsonCache
from .tle import get_TLE, Tle
from .core import predict_single_satellite_overpasses
from ._time import julian_date
from .satellites import SatellitePredictor
from .locations import Location
from .sources import TLE


@click.command()
@click.option('-s', '--satellite-id', type=int)  # satellite id
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
def main(satellite_id, utc_offset, days, latitude, longitude, height, twelve, alltypes, quiet, verbose, no_cache):
    """
    Command line interface for pass predictions
    """
    tz = datetime.timezone(datetime.timedelta(hours=utc_offset))
    location = Location(
        latitude_deg=latitude,
        longitude_deg=longitude,
        elevation_m=height,
        name=""
    )
    satellite = Satellite(
        id=satellite_id,
    )
    date_start = datetime.datetime.utcnow()
    min_elevation = 10.0 # degrees

    if no_cache:
        tle = get_TLE(satellite.id)
    else:
        with JsonCache() as cache:
            sat_key = f'sat:{satellite.id}'
            res = cache.get(sat_key)
            if res:
                tle = Tle(tle1=res['tle1'], tle2=res['tle2'])
            else:
                tle = get_TLE(satellite.id)
                cache.set(sat_key, tle.dict(), ttl=86400)

    predictor = SatellitePredictor.from_tle(tle)
    overpasses = predict_single_satellite_overpasses(predictor, location, date_start, days, min_elevation=min_elevation)
    overpass_table(overpasses, location, tle, tz, twelvehour=twelve, quiet=quiet)
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
    table_header +=  f"            {'Start':^17s}   {'Maximum':^17s}   {'End':^17s}\n"
    # table_header +=  "  Date    Mag   {0}   {0}   {0}      Type\n".format(point_header)
    table_header +=  "  Date    {0}   {0}   {0}\n".format(point_header)
    table_header += "--------  "
    table_header += point_header_underline + " "*3
    table_header += point_header_underline + " "*3
    table_header += point_header_underline + " "*3
    # table_header += "-"*10   # for 'Type' underline
    click.echo(table_header)

    def point_string(point):
        time = point.dt.astimezone(tz)
        if twelvehour:
            point_line = time.strftime("%I:%M:%S") + time.strftime("%p")[0].lower()
        else:
            point_line = time.strftime("%H:%M:%S")
        point_line += " " + "{:>2}\u00B0".format(int(point.elevation))
        point_line += " " + "{:3}".format(point.direction)
        return point_line

    for overpass in overpasses:
        table_data = "{}".format(overpass.aos.dt.astimezone(tz).strftime("%m/%d/%y"))
        # if overpass.brightness is not None:
        #     brightness_str = f"{overpass.brightness:4.1f}"
        # else:
        #     brightness_str = " "*4
        # table_data += " "*2 + brightness_str
        table_data += " "*2 + point_string(overpass.aos) + ' |'
        table_data += " " + point_string(overpass.tca) + ' |'
        table_data += " " + point_string(overpass.los)
        if overpass.type:
            table_data += " "*4 + f'{overpass.type.value:^9}'# + '\n'
            fg = 'green' if overpass.type.value == PassType.visible else None
        else:
            fg = None
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