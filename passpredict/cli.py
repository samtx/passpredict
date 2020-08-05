import datetime
from pprint import pprint

import click

from .schemas import Satellite, Tle, Location, PassType
from .timefn import truncate_datetime, get_date
from .geocoding import geocoder
from .utils import get_TLE, Cache
from .predictions import predict

@click.command()
@click.option('-s', '--satellite-id', type=int, prompt=True)  # satellite id
@click.option('-l', '--location', 'location_string', type=str, prompt=True)  # location string
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
    # Prompt for location
    # query = input('Enter location: ')
    data = geocoder(location_string)

    # tz_offset_str = input('Enter timezone offset: ')
    # tz_offset = float(tz_offset_str)
    tz = datetime.timezone(datetime.timedelta(hours=utc_offset))

    location = Location(
        lat=float(data['lat']),
        lon=float(data['lon']),
        h=0.0,
        name=location_string
    )
    satellite = Satellite(
        id=satellite_id,
        # name="Int. Space Station"
    )
    tle = get_TLE(satellite)
    date_start = datetime.date.today()
    date_end = date_start + datetime.timedelta(days=days)
    min_elevation = 10.01 # degrees
    cache = Cache() if not no_cache else None
    overpasses = predict(location, satellite, date_start=date_start, date_end=date_end, dt_seconds=1, min_elevation=min_elevation, verbose=verbose, cache=cache)
    table_str = overpass_table(overpasses, location, tle, tz, twelvehour=twelve, alltypes=alltypes, quiet=quiet)
    print(table_str)
    return 0


def overpass_table(overpasses, location, tle, tz=None, twelvehour=False, alltypes=False, quiet=False):
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
    table_data = ""
    if not quiet:
        satellite_id = tle.satellite.id
        # Print datetimes with the correct timezone
        table_title += f"Satellite ID {satellite_id} overpasses for {location.name:s}\n"
        table_title += f"Lat={location.lat:.4f}\u00B0, Lon={location.lon:.4f}\u00B0, Timezone {tz}\n"
        table_title += f"Using TLE\n"
        table_title += f"{tle.tle1:s}\n"
        table_title += f"{tle.tle2:s}\n\n"
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
    table_header += "\n"

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
        if not alltypes:
            # filter overpasses to visible only
            if overpass.type != PassType.visible:
                continue
        table_data += "{}".format(overpass.start_pt.datetime.astimezone(tz).strftime("%m/%d/%y"))
        if overpass.brightness is not None:
            brightness_str = f"{overpass.brightness:4.1f}"
        else:
            brightness_str = " "*4
        table_data += " "*2 + brightness_str
        table_data += " "*2 + point_string(overpass.start_pt) + ' |'
        table_data += " " + point_string(overpass.max_pt) + ' |'
        table_data += " " + point_string(overpass.end_pt)
        table_data += " "*4 + f'{overpass.type.value:^9}'
        table_data += "\n"
    return table_title + table_header + table_data
# --------- --- ---

def plot_elevation(date, elevation):
    import matplotlib.pyplot as plt
    plt.plot(date, elevation)
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()