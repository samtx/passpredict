import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
import datetime
from passpredict.rotations import rot1, rot2, rot3, theta_GMST1982, site_sat_rotations
from passpredict.solar import sun_pos, is_sat_illuminated
from passpredict.topocentric import razel
from passpredict.constants import (
    R_EARTH, R2_EARTH, e_EARTH, e2_EARTH, MU, J2, J2000, AU_M, AU_KM, ASEC360,
    DAY_S, ASEC2RAD, DEG2RAD, RAD2DEG, tau
)
from passpredict.propagate import propagate, get_TLE
from passpredict.timefn import jday2datetime
from passpredict.models import Point, Overpass
import pickle

def vector_angle(r1, r2):
    """Compute the angle between two vectors

        r1 and r2 : (3, n), n is the number of observations
    """
    numerator = np.einsum("ij,ij->j", r1, r2)  # vectorized dot product
    denominator = norm(r1, axis=0) * norm(r2, axis=0)
    out = np.arccos(numerator / denominator) * RAD2DEG
    return out


def satellite_visible(rsatECI, rsiteECI, rho, jdt):
    """Determine visibility of satellite from Earth"""
    visible = np.zeros(jdt.size)
    # First, determine if satellite is above the horizon
    # find indecies where rho[2] > 0 --> idx
    vis_idx = np.nonzero(rho[2] > 0)[0]
    # Loop over indecies for times that satellite is over horizon
    for i in range(len(vis_idx)):
        idx = vis_idx[i]
        jdt_i = jdt[idx]
        rsatECI_i = rsatECI[:, idx]
        rsiteECI_i = rsiteECI[:, idx]
        rho_i = rho[:, idx]
        # Check if site is in daylight, compute dot product of sun position.
        rsun = sun_pos(jdt_i)
        if len(rsun.shape) > 1:
            rsun = rsun.flatten()
        # TODO: compute ECI vector for site position
        if np.dot(rsun, rsiteECI) > 0:
            # site is in daylight
            visible[idx] = 1
        else:
            # If nighttime, check if satellite is illuminated or in shadow
            if is_sat_illuminated(rsatECI_i, rsun):
                # satellite is illuminated
                visible[idx] = 3
            else:
                # satellite is in shadow
                visible[idx] = 2
    return visible


def get_overpasses(el, azm, rng, jdt_ary, rSEZ, rsiteECI=None, rsatECI=None, min_elevation=10, loc=None, sat=None):
    # # change julian dates to datetimes
    # num_jdt = jdt_ary.size
    # dt_array = np.empty(num_jdt, dtype=object)
    # for i in range(num_jdt):
    #     dt_array[i] = jday2datetime(jdt_ary[i])
    el0 = el[:-1] - min_elevation
    el1 = el[1:] - min_elevation
    el_change_sign = (el0*el1 < 0)
    # Find the start of an overpass
    start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]
    # Find the end of an overpass
    end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]
    # Iterate over start/end indecies and gather inbetween indecies
    num_overpasses = min(start_idx.size, end_idx.size)
    if start_idx.size < end_idx.size:
        end_idx = end_idx[1:]
    overpasses = np.empty(num_overpasses, dtype=object)
    for j in range(num_overpasses):
        # Store indecies of overpasses in a list
        idx0 = start_idx[j]
        idxf = end_idx[j]
        overpass_idx = np.arange(idx0, idxf+1, dtype=int)
        idxmax = np.argmax(el[overpass_idx])

        start_pt = Point(
            jday2datetime(jdt_ary[idx0]),
            azm[idx0],
            el[idx0],
            rng[idx0]
        )
        max_pt = Point(
            jday2datetime(jdt_ary[idx0 + idxmax]),
            azm[idx0 + idxmax],
            el[idx0 + idxmax],
            rng[idx0 + idxmax]
        )
        end_pt = Point(
            jday2datetime(jdt_ary[idxf]),
            azm[idxf],
            el[idxf],
            rng[idxf]
        )
        # sat_vis = satellite_visible(rsatECI, rsiteECI, rSEZ, jdt)
        overpass = Overpass(
            loc,
            sat,
            start_pt,
            max_pt,
            end_pt,
            jdt_ary,
            rSEZ[:,overpass_idx]
        )
        overpasses[j] = overpass
    return overpasses


def predict_passes(lat, lon, h, rsatECEF, rsatECI, jdt, rsun=None, min_elevation=None, loc=None, sat=None):
    rSEZ = site_sat_rotations(lat, lon, h, rsatECEF)
    # rsiteECI = site2eci(lat, lon, h, jdt)
    rng, az, el = razel(rSEZ)
    # plot_elevation(np.arange(el.size), el)
    overpasses = get_overpasses(el, az, rng, jdt, rSEZ, rsiteECI=None, rsatECI=None, min_elevation=min_elevation, loc=loc, sat=sat)
    return overpasses


def predict(location, satellite, dt_start=None, dt_end=None, dt_seconds=1, min_elevation=None, tle=None, reload=True):
    """
    Full prediction algorithm:
      1. Download TLE data
      2. Propagate satellite using SGP4
      3. Predict overpasses based on site location
      4. Return overpass object and print to screen

    Params:
        location : Location object
            latitude of site location, in decimal, north is positive
        satellite: Satellite object
            satellite ID number in Celestrak, ISS is 25544
    """
    if dt_start is None:
        dt_start = datetime.datetime.now()
    if dt_end is None:
        dt_end = dt_start + datetime.timedelta(days=14)
    if tle is None:
        tle = get_TLE(satellite)
    print(f"begin propagation from {dt_start.isoformat()} to {dt_end.isoformat()}")
    _reload = reload
    if _reload:
        satellite_rv = propagate(tle.tle1, tle.tle2, dt_start, dt_end, dt_seconds)
        satellite_rv.satellite = satellite
        satellite_rv.tle = tle
        with open(f'satellite_{satellite.id:d}.pkl', 'wb') as f:
            pickle.dump(satellite_rv, f)
    else:
        with open(f'satellite_{satellite.id:d}.pkl', 'rb') as f:
            satellite_rv = pickle.load(f)
    print('begin prediction...')
    # set minimum elevation parameter: min_elevation = 10 degrees
    overpasses = predict_passes(
        location.lat, location.lon, location.h,
        satellite_rv.rECEF, satellite_rv.rECI, satellite_rv.julian_date,
        min_elevation=min_elevation)
    return overpasses


def overpass_table(overpasses, location, tle, twentyfourhour=True):
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
    tz = location.tz
    # Print datetimes with the correct timezone
    table_title = ""
    table_title += f"{satellite_name:s} overpasses for {location.name:s}\n"
    table_title += f"Lat={location.lat:.4f}\u00B0, Lon={location.lon:.4f}\u00B0, Timezone {tz.zone}\n"
    table_title += f"Using TLE\n"
    table_title += f"{tle.tle1:s}\n"
    table_title += f"{tle.tle2:s}\n\n"
    if not twentyfourhour:
        point_header = "  Time    El\u00B0 Az\u00B0"
        point_header_underline = "--------- --- ---"
    else:
        point_header = "  Time   El\u00B0 Az\u00B0"
        point_header_underline = "-------- --- ---"
    table_header =  "Date     St[{0}] Mx[{0}] En[{0}]\n".format(point_header)
    table_header += "--------    "
    table_header += point_header_underline + " "*5
    table_header += point_header_underline + " "*5
    table_header += point_header_underline + " "*5
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
        table_data += " "*4 + point_string(overpass.start_pt)
        table_data += " "*5 + point_string(overpass.max_pt)
        table_data += " "*5 + point_string(overpass.end_pt)
        table_data += "\n"
    return table_title + table_header + table_data


def plot_elevation(date, elevation):
    import matplotlib.pyplot as plt
    plt.plot(date, elevation)
    plt.grid()
    plt.show()


if __name__ == "__main__":
    from passpredict.models import Location, Satellite, Tle
    from passpredict.timefn import truncate_datetime
    from pprint import pprint
    import pytz

    # tle1 = "1 25544U 98067A   20144.05946705  .00016717  00000-0  10270-3 0  9060"
    # tle2 = "2 25544  51.6397 112.8409 0001237 325.9375  34.1696 15.49401715 28113"
    tle1 = "1 25544U 98067A   20145.02695725  .00016717  00000-0  10270-3 0  9072"
    tle2 = "2 25544  51.6412 108.0531 0001249 337.3712  22.7382 15.49390188 28261"
    epoch = datetime.datetime(2020, 5, 23, 1, 25, 37, tzinfo=pytz.utc)  # 23 May 2020 01:25:37 UTC
    austin = Location(30.2711, -97.7437, 0.0, 'Austin', pytz.timezone('US/Central'))
    satellite = Satellite(25544, "Int. Space Station")
    tle = Tle(tle1, tle2, epoch, satellite)
    dt_start = truncate_datetime(datetime.datetime(2020,5,21,0,0,0))# - datetime.timedelta(days=1)
    dt_end = dt_start + datetime.timedelta(days=5)
    min_elevation = 10.1 # degrees
    overpasses = predict(austin, satellite, dt_start=dt_start, dt_end=dt_end, dt_seconds=1, min_elevation=min_elevation, reload=True)
    print('begin printing table...')
    table_str = overpass_table(overpasses, austin, tle, False)
    print(table_str)

