from datetime import datetime, timedelta, timezone
from functools import cached_property
import json
from typing import NamedTuple, Float

import requests
import numpy as np
from pydantic import BaseModel, Field

from .utils import grouper


class TleSchema(BaseModel):
    tle1: str
    tle2: str
    epoch: datetime
    satid: int

    class Config:
        title = 'TLE'


class Tle():
    def __init__(self, tle1: str, tle2: str):
        self.tle1 = tle1
        self.tle2 = tle2
        
    @cached_property
    def epoch(self) -> datetime:
        return epoch_from_tle(self.tle1)

    @cached_property
    def satid(self) -> int:
        return satid_from_tle(self.tle1)

    def to_schema(self):
        return TleSchema(
            tle1=self.tle1,
            tle2=self.tle2,
            epoch=self.epoch,
            satid=self.satid    
        )


def epoch_from_tle_datetime(epoch_string: str) -> datetime:
    """
    Return datetime object from tle epoch string
    """
    epoch_year = int(epoch_string[0:2])
    if epoch_year < 57:
        epoch_year += 2000
    else:
        epoch_year += 1900
    epoch_day = float(epoch_string[2:])
    epoch_day, epoch_day_fraction = np.divmod(epoch_day, 1)
    epoch_microseconds = epoch_day_fraction * 24 * 60 * 60 * 1e6
    return datetime(epoch_year, month=1, day=1, tzinfo=timezone.utc) + \
        timedelta(days=int(epoch_day-1)) + \
        timedelta(microseconds=int(epoch_microseconds)
    )
    

def epoch_from_tle(tle1: str) -> datetime:
    """
    Extract epoch as datetime from tle line 1
    """
    epoch_string = tle1[18:32]
    return epoch_from_tle_datetime(epoch_string)

def epoch_to_jd(tle1: str) -> float:
    """
    Convert TLE epoch to julian date
    """
    epoch_string = tle1[18:32]
    yr = int(epoch_string[0:2])
    if yr < 57:
        yr += 2000
    else:
        yr += 1900
    dy = float(epoch_string[2:])
    dy, fr = divmod(dy, 1)
    us = fr * 24 * 60 * 60 * 1e6
    # def julian_date(yr, mo=None, dy=None, hr=None, mn=None, sec=None, us=0.0):
    yr, mo, dy = dt.year, dt.month, dt.day
    hr, mn, sec = dt.hour, dt.minute, dt.second
    sec += dt.microsecond * (10 ** -6)
    sec += us * 10e-6
    jd1 = 367 * yr
    jd2 = 7 * (yr + (mo + 9) // 12) // 4
    jd3 = (275 * mo) // 9
    jd4 = dy
    jd5 = 1721013.5
    jd6 = ((sec / 60 + mn) / 60 + hr) / 24
    jd = jd1 - jd2 + jd3 + jd4 + jd5 + jd6
    # print([jd1, jd2, jd3, jd4, jd5, jd6])
    return jd

    

def satid_from_tle(tle1: str) -> int:
    """
    Extract satellite NORAD ID as int from tle line 1
    """
    return int(tle1[2:7])


def get_orbit_data_from_celestrak(satellite_id):
    """

    Params:
        satellite_id : int
            NORAD satellite ID


    See https://celestrak.com/NORAD/documentation/gp-data-formats.php

    Can use the new celestrak api for satellite ID
    https://celestrak.com/NORAD/elements/gp.php?CATNR=25544&FORMAT=json

    for tle api:
    https://celestrak.com/satcat/tle.php?CATNR=25544

    Supplemental TLEs available: (not fully working as json)
    https://celestrak.com/NORAD/elements/supplemental/gp-index.php?GROUP=iss&FORMAT=json

    https://celestrak.com/NORAD/elements/supplemental/starlink.txt
    https://celestrak.com/NORAD/elements/supplemental/iss.txt
    
    """
    query = {
        'CATNR': satellite_id,
        'FORMAT': 'json'
    }
    url = 'https://celestrak.com/NORAD/elements/gp.php'
    r = requests.get(url, data=query)
    return r.json()


def parse_tles_from_celestrak(satellite_id=None):
    """
    Download current TLEs from Celestrak and save them to a JSON file
    
    """
    if satellite_id is None:
        url = 'https://celestrak.com/NORAD/elements/stations.txt'
        params = {}
    else:
        url = 'https://celestrak.com/satcat/tle.php'
        params = {'CATNR': satellite_id}
    r = requests.get(url, params=params, stream=True)
    tle_data = {}
    for tle_strings in grouper(r.text.splitlines(), 3):
        tle_data.update(parse_tle(tle_strings))
    return tle_data


def parse_tle(tle_string_list):
    """
    Parse a single 3-line TLE from celestrak
    """
    tle0, tle1, tle2 = tle_string_list
    name = tle0.strip()  # satellite name
    satellite_id = satid_from_tle(tle1)
    return {satellite_id : {'name': name, 'tle1': tle1, 'tle2': tle2}}


def get_TLE(satid: int, tle_data=None):
    tle_data = parse_tles_from_celestrak(satid)
    tle1 = tle_data[satid]['tle1']
    tle2 = tle_data[satid]['tle2']
    return Tle(tle1=tle1, tle2=tle2)
    

def save_TLE_data(url=None):
    tle_data = parse_tles_from_celestrak(url)
    with open('tle_data.json', 'w') as file:
        json.dump(tle_data, file)

class OMM(NamedTuple):
    """
    Data structure for holding orbital elements
    """
    jdsatepoch: float   # julian date
    jdsatepochF: float  # julian date fraction
    no_kozai: float     # kozai mean motion [rev/day], line 2, ch 53-63
    ecco: float         # eccentricity, line 2, ch 27-33
    inclo: float        # inclination [deg], line 2, ch 9-16
    nodeo: float        # right ascension of ascending node [deg], line 2, ch 18-25
    argpo: float        # argument of perigee [deg] line 2, ch 35-42
    mo: float           # mean anomolay, line 2, ch 44-51
    nddot: float        # second derivative of mean motion, line 1
    bstar: float        # B star drag term, line 1
    ndot: float         # first derivative of mean motion, line 1
    elnum: int          # element number, line 1
    revnum: int         # revolution number at epoch
    classification: str = 'U'
    ephtype: int = 0    # element set type


def tle_to_omm(tle1: str, tle2: str) -> OMM:
    """
    Convert TLE strings to OMM data 
    """


#     // sgp4fix demonstrate method of running SGP4 directly from orbital element values
# 	//1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
# 	//2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656
# 	const double deg2rad = pi / 180.0;         //   0.0174532925199433
# 	const double xpdotp = 1440.0 / (2.0 *pi);  // 229.1831180523293

# 	whichconst = wgs72;
# 	opsmode = 'a';
# 	// new alpha5 or 9-digit number
# #ifdef _MSC_VER
# 	strcpy_s(satrec.satnum, sizeof(satrec.satnum), "8195");
# #else
# 	strcpy(satrec.satnum, "8195");
# #endif

# 	satrec.jdsatepoch = 2453911.0;
# 	satrec.jdsatepochF = 0.8321544402;
# 	satrec.no_kozai = 2.00491383;
# 	satrec.ecco = 0.6877146;
# 	satrec.inclo = 64.1586;
# 	satrec.nodeo = 279.0717;
# 	satrec.argpo = 264.7651;
# 	satrec.mo = 20.2257;
# 	satrec.nddot = 0.00000e0;
# 	satrec.bstar = 0.11873e-3;
# 	satrec.ndot = 0.00000099;
# 	satrec.elnum = 813;
# 	satrec.revnum = 22565;
# 	satrec.classification = 'U';
# 	strncpy_s(satrec.intldesg, "          ", 11 * sizeof(char));
# 	satrec.ephtype = 0;

# 	// convert units and initialize
# 	satrec.no_kozai = satrec.no_kozai / xpdotp; //* rad/min
# 	satrec.ndot = satrec.ndot / (xpdotp*1440.0);  //* ? * minperday
# 	satrec.nddot = satrec.nddot / (xpdotp*1440.0 * 1440);
# 	satrec.inclo = satrec.inclo  * deg2rad;
# 	satrec.nodeo = satrec.nodeo  * deg2rad;
# 	satrec.argpo = satrec.argpo  * deg2rad;
# 	satrec.mo = satrec.mo     * deg2rad;

# 	// set start/stop times for propagation
# 	startmfe = 0.0;
# 	stopmfe = 2880.0;
# 	deltamin = 120.0;

# 	SGP4Funcs::sgp4init(whichconst, opsmode, satrec.satnum, satrec.jdsatepoch + satrec.jdsatepochF - 2433281.5, satrec.bstar,
# 		satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
# 		satrec.nodeo, satrec);