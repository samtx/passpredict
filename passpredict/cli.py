import datetime
from pprint import pprint

import click

from .schemas import Satellite, Tle, Location
from .timefn import truncate_datetime
from .geocoding import geocoder
from .utils import get_TLE
from .predictions import predict

@click.command()
@click.option('-12/-24', default=False)  # 12 hour / 24 hour format
@click.option('-u', '--utc-offset', type=float)  # utc offset
@click.option('-d', '--days', type.click.IntRange(1, 14, clamp=True)) # day range
def main():
    """
    Command line interface for pass predictions
    """
    # Prompt for location
    query = input('Enter location: ')
    data = geocoder(query)

    tz_offset_str = input('Enter timezone offset: ')
    tz_offset = float(tz_offset_str)
    tz = datetime.timezone(tz_offset)

    location = Location(
        lat=float(data['lat']),
        lon=float(data['lon']),
        h=0.0,
        name=query
    )
    satellite = Satellite(
        id=25544,
        name="Int. Space Station"
    )
    tle = get_TLE(satellite)
    dt_start = truncate_datetime(datetime.datetime.now())# - datetime.timedelta(days=1)
    dt_end = dt_start + datetime.timedelta(days=10)
    min_elevation = 10.01 # degrees
    overpasses = predict(location, satellite, dt_start=dt_start, dt_end=dt_end, dt_seconds=5, min_elevation=min_elevation, reload=True)
    print('begin printing table...')
    table_str = overpass_table(overpasses, location, tle, tz, twentyfourhour=True)
    print(table_str)

def overpass_table(overpasses, location, tle, tz=None, twentyfourhour=True):
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
    satellite_name = tle.satellite.name
    # Print datetimes with the correct timezone
    table_title = ""
    table_title += f"{satellite_name:s} overpasses for {location.name:s}\n"
    table_title += f"Lat={location.lat:.4f}\u00B0, Lon={location.lon:.4f}\u00B0, Timezone {tz}\n"
    table_title += f"Using TLE\n"
    table_title += f"{tle.tle1:s}\n"
    table_title += f"{tle.tle2:s}\n\n"
    if not twentyfourhour:
        point_header = "  Time    El\u00B0 Az\u00B0"
        point_header_underline = "--------- --- ---"
    else:

        #   Time   Elx Azx
        point_header = "  Time   El\u00B0 Az\u00B0"
        point_header_underline = "-------- --- ---"
    table_header =  f"           {'Start':^17s}   {'Maximum':^17s}   {'End':^17s}\n"
    table_header +=  "  Date     {0}   {0}   {0}\n".format(point_header)
    table_header += "--------   "
    table_header += point_header_underline + " "*3
    table_header += point_header_underline + " "*3
    table_header += point_header_underline + " "*3
    table_header += "\n"

    def point_string(point):
        time = point.datetime.astimezone(tz)
        if not twentyfourhour:
            point_line = time.strftime("%I:%M:%S") + time.strftime("%p")[0].lower()
        else:
            point_line = time.strftime("%H:%M:%S")
        point_line += " " + "{:>3}".format(int(point.elevation))
        point_line += " " + "{:3}".format(point.direction_from_azimuth())
        return point_line

    table_data = ""
    for overpass in overpasses:
        table_data += "{}".format(overpass.start_pt.datetime.astimezone(tz).strftime("%m/%d/%y"))
        table_data += " "*3 + point_string(overpass.start_pt) + ' |'
        table_data += " " + point_string(overpass.max_pt) + ' |'
        table_data += " " + point_string(overpass.end_pt)
        table_data += "\n"
    return table_title + table_header + table_data
# --------- --- ---

def plot_elevation(date, elevation):
    import matplotlib.pyplot as plt
    plt.plot(date, elevation)
    plt.grid()
    plt.show()
