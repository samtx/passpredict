import json
from datetime import datetime, timedelta, timezone, tzinfo
from pathlib import Path
import collections
import itertools
import typing

import tqdm

import astrodynamics as astro


class JSONEncoder(json.JSONEncoder):
    """
    Custom json encoder for seializing datetime objects in isoformat
    """
    def default(self, obj):
        if isinstance(obj, datetime.datetime):
            return obj.isoformat()
        return super().default(obj)


def main():
    """
    Generate pass predictions for the next 4 days for the given locations and satellites
    """
    Location = collections.namedtuple('Location', 'name, lat, lon, h')
    Satellite = collections.namedtuple('Satellite', 'name, satid, lines')
    locations = [
        Location('Svalbard Satellite Station, Norway', 78.2297, 15.3975, 458),
        Location('Inuvik, Canada', 68.3607, -133.7230, 15),
        Location('Kiev, Ukraine', 50.4501, 30.5234, 179),
        Location('Austin, Texas, USA', 32.1234, -97.1234, 0),
        Location('Quito, Ecuador', 0.1807, -78.4678, 2850),
        Location('Johannesburg, South Africa', -26.2041, 28.0473, 1753),
        Location('New Norcia Deep Space Station, Australia', -31.0483, 116.1914, 252),
        Location('McMurdo Station, Antarctica', -77.8419, 166.6863, 10),
    ]
    satellites = [
        Satellite('International Space Station', 25544, ('1 25544U 98067A   21274.89460679  .00005555  00000-0  10931-3 0  9992','2 25544  51.6449 175.1112 0004190  48.8354  53.9444 15.48895782305133')),
        Satellite('Hubble Space Telescope', 20580, ('1 20580U 90037B   21274.54126351  .00001079  00000-0  54194-4 0  9991','2 20580  28.4698  83.0715 0002784 122.8787 269.0686 15.09757146527044')),
        Satellite('Terra', 25994, ('1 25994U 99068A   21274.75027490  .00000107  00000-0  33588-4 0  9990','2 25994  98.1678 347.1838 0001354  91.6968  61.5062 14.57132833158940')),
        Satellite('Envisat', 27386, ('1 27386U 02009A   21274.71722274  .00000080  00000-0  39878-4 0  9992','2 27386  98.1579 263.3896 0001195  87.0592  86.1444 14.38044645 26342')),
        Satellite('Intelsat 5', 24916, ('1 24916U 97046A   21274.50960067  .00000091  00000-0  00000-0 0  9997','2 24916   6.6785  62.7523 0003074 141.8639 212.2992  1.00271606 88373')),
        Satellite('Aprizesat 5', 37792, ('1 37792U 11044E   21274.60878477  .00000190  00000-0  35454-4 0  9996','2 37792  98.3609 136.8622 0059209 108.1425 252.6243 14.75523908544976')),
        Satellite('OFEQ 11', 41759, ('1 41759U 16056A   20133.81295366 0.00000000  00000-0  00000-0 0    06','2 41759 141.8552 279.8922 0006334 285.1363  74.8629 15.34726396    09'))
    ]
    start = datetime(2021, 10, 2, 0, 0, 0, tzinfo=timezone.utc)
    end = datetime(2021, 10, 6, 0, 0, 0, tzinfo=timezone.utc)
    data = {
        'locations': [l._asdict() for l in locations],
        'satellites': [s._asdict() for s in satellites],
        'start': start.isoformat(),
        'end': end.isoformat(),
        'overpasses': {},
    }
    loc_sat = itertools.product(locations, satellites)


    fname = Path(__file__).parent / 'predictions.json'
    with open(fname, 'w') as f:
        json.dump(data, f, indent=2, cls=JSONEncoder)



def find_overpasses(location, satellite) -> typing.List[astro.PredictedPass]:
    """
    Find overpasses
    """
    pass

if __name__ == "__main__":
    main()