import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
import datetime
from passpredictor.rotations import rot1, rot2, rot3, theta_GMST1982, site_sat_rotations
from passpredictor.solar import sun_pos, is_sat_illuminated
from passpredictor.topocentric import razel
from passpredictor.constants import (
    R_EARTH, R2_EARTH, e_EARTH, e2_EARTH, MU, J2, J2000, AU_M, AU_KM, ASEC360,
    DAY_S, ASEC2RAD, DEG2RAD, RAD2DEG, tau
)
from passpredictor.propagate import propagate, get_TLE
from passpredictor.timefn import jday2datetime
from passpredictor.models import Point, Overpass


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


def get_overpasses(el, azm, rng, jdt_ary, rSEZ, rsiteECI=None, rsatECI=None, loc=None, sat=None):
    # # change julian dates to datetimes
    # num_jdt = jdt_ary.size
    # dt_array = np.empty(num_jdt, dtype=object)
    # for i in range(num_jdt):
    #     dt_array[i] = jday2datetime(jdt_ary[i])

    el0 = el[:-1]
    el1 = el[1:]
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


def predict_passes(lat, lon, h, rsatECEF, rsatECI, jdt, rsun=None, loc=None, sat=None):
    rSEZ = site_sat_rotations(lat, lon, h, rsatECEF)
    # rsiteECI = site2eci(lat, lon, h, jdt)
    rng, az, el = razel(rSEZ)
    plot_elevation(np.arange(el.size), el)
    overpasses = get_overpasses(el, az, rng, jdt, rSEZ, rsiteECI=None, rsatECI=None, loc=loc, sat=sat)
    return overpasses


def predict(lat, lon, h, satellite, dt_begin=None, dt_end=None, dt_seconds=1):
    """
    Full prediction algorithm:
      1. Download TLE data
      2. Propagate satellite using SGP4
      3. Predict overpasses based on site location
      4. Return overpass object and print to screen

    Params:
        lat : float
            latitude of site location, in decimal, north is positive
        lon : float
            longitude of site location, in decimal, east is positive
        h : float
            elevation of site in meters
        satid: int
            satellite ID number in Celestrak, ISS is 25544

    """
    if dt_begin is None:
        dt_begin = datetime.datetime.now()
    if dt_end is None:
        dt_end = dt_begin + datetime.timedelta(days=14)
    tle = get_TLE(satellite)
    print(f'begin propagation from {dt_begin} to {dt_end}')
    satellite_rv = propagate(tle.tle1, tle.tle2, dt_begin, dt_end, dt_seconds)
    satellite_rv.satellite = satellite
    satellite_rv.tle = tle
    print('begin prediction...')
    overpasses = predict_passes(lat, lon, h, satellite_rv.rECEF, satellite_rv.rECI, satellite_rv.julian_date)
    return overpasses


def overpass_table(overpasses=None):
    """
    Return a formatted string for tabular output

    Params:
        overpasses: list
            A list of Overpass objects

    Return:
        table : str
            tabular formatted string
    """
    point_header = "  Time     Az\u00B0   El\u00B0"
    point_header_underline = "--------  ----  ----"
    header =  "Date      St[{0}]  Mx[{0}]  En[{0}]\n".format(point_header)
    header_length = len(header)
    header += "--------     "
    header += point_header_underline + " "*6
    header += point_header_underline + " "*6
    header += point_header_underline + " "*6
    print(header)
    # print(header_length)

    def point_string(point):
        point_line = point.datetime.strftime("%I:%M:%S")
        point_line += " "*2 + "{:>4}".format(int(point.azimuth))
        point_line += " "*2 + "{:>4}".format(int(point.elevation))
        return point_line

    for overpass in overpasses:
        line = "{}".format(overpass.start_pt.datetime.strftime("%m/%d/%y"))
        line += " "*5 + point_string(overpass.start_pt)
        line += " "*6 + point_string(overpass.max_pt)
        line += " "*6 + point_string(overpass.end_pt)
        print(line)

def plot_elevation(date, elevation):
    import matplotlib.pyplot as plt
    plt.plot(date, elevation)
    plt.grid()
    plt.show()


if __name__ == "__main__":


    from passpredictor.models import Location, Satellite
    from pprint import pprint
    austin = Location(30.2672, -97.7431, 0.0, 'Austin')
    satellite = Satellite(25544, "Int. Space Station")
    dt_end = datetime.datetime.now() + datetime.timedelta(days=1)
    overpasses = predict(austin.lat, austin.lon, austin.h, satellite, dt_end=dt_end)
    print('begin printing table...')
    overpass_table(overpasses)

