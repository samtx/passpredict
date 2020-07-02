import json
import datetime
import os
from itertools import zip_longest
from collections import OrderedDict

import numpy as np
import requests

from .schemas import Tle

CACHE_DIRECTORY = ".passpredict_cache"


def shift_angle(x: float) -> float:
    """Shift angle in radians to [-pi, pi)
    
    Args:
        x: float, angle in radians

    Reference: 
        https://stackoverflow.com/questions/15927755/opposite-of-numpy-unwrap/32266181#32266181
    """
    return (x + np.pi) % (2 * np.pi) - np.pi
    


def grouper(iterable, n, fillvalue=None):
    """
    from itertools recipes https://docs.python.org/3.7/library/itertools.html#itertools-recipes
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def epoch_from_tle_datetime(epoch_string: str) -> datetime.datetime:
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
    epoch = datetime.datetime(epoch_year, month=1, day=1) + \
            datetime.timedelta(days=int(epoch_day-1)) + \
            datetime.timedelta(microseconds=int(epoch_microseconds))
    return epoch
    

def epoch_from_tle(tle1: str) -> datetime.datetime:
    """
    Extract epoch as datetime from tle line 1
    """
    epoch_string = tle1[18:32]
    return epoch_from_tle_datetime(epoch_string)
    

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
    satellite_id = tle1[2:7]
    return {satellite_id : {'name': name, 'tle1': tle1, 'tle2': tle2}}


def get_TLE(satellite, tle_data=None):
    if tle_data is None:
        if not os.path.exists('tle_data.json'):
            tle_data = parse_tles_from_celestrak()
            with open('tle_data.json', 'w') as file:
                json.dump(tle_data, file)
        else:
            with open('tle_data.json', 'r') as file:
                tle_data = json.load(file)
    tle1 = tle_data[str(satellite.id)]['tle1']
    tle2 = tle_data[str(satellite.id)]['tle2']
    epoch = epoch_from_tle(tle1)
    tle = Tle(tle1=tle1, tle2=tle2, epoch=epoch, satellite=satellite)
    return tle


def save_TLE_data(url=None):
    tle_data = parse_tles_from_celestrak(url)
    with open('tle_data.json', 'w') as file:
        json.dump(tle_data, file)


class Cache():
    def __init__(self, cache_directory=CACHE_DIRECTORY, ttl=30):
        self.directory = cache_directory
        self.ttl = ttl # days
        self.cache = OrderedDict()

    def set(self, key, val):
        pass

    def get(self, key, val):
        return None


def cache_tle(data, method, cache_directory=CACHE_DIRECTORY):
    """
    Get/Set satellite TLE data

    Params:
        data = tle data
        method = str ['get','set']
        cache_directory: str, optional
    """
    pass


def cache_satellite_position(data, method, cache_directory=CACHE_DIRECTORY):
    """
    Get/Set satellite ECEF position vectors in cache
    """
    pass
