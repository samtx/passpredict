import math
from datetime import datetime, timedelta, timezone
from functools import cached_property
import json
from typing import NamedTuple, Union, Tuple
import dataclasses
from itertools import zip_longest

import httpx
import numpy as np

from passpredict._time import epoch_to_jd, jday2datetime_us


# class TleSchema(BaseModel):
#     tle1: str
#     tle2: str
#     epoch: datetime
#     satid: int

#     class Config:
#         title = 'TLE'


# from orbit_predictor.sources
class TLE(NamedTuple):
    satid: Union[int, str]        # NORAD satellite ID
    lines: Tuple[str]   # tuple of tle strings (tle1, tle2)
    name: str = ""

    @property
    def sate_id(self):
        return self.satid

    @property
    def epoch(self) -> datetime:
        return epoch_from_tle(self.lines[0])

    @property
    def tle1(self) -> str:
        return self.lines[0]

    @property
    def tle2(self) -> str:
        return self.lines[1]

    def dict(self):
        return self._asdict()

    def __repr__(self):
        return f"<TLE satid={self.satid} epoch={self.epoch}>"


def jd_to_epoch_string(jd: float) -> str:
    """
    Get TLE epoch string from julian date
    """
    d = jday2datetime_us(jd)
    yearday = d.strftime("%y%j")
    dayfraction = (d.hour*3600 + d.minute*60 + d.second + d.microsecond*1e-6) / 86400
    dayfraction = round(dayfraction * 1e8, 0)
    s = f"{yearday}.{dayfraction:8.0f}"
    return s


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


# class OMM(NamedTuple):
#     """
#     Data structure for holding orbital elements
#     """
#     jdsatepoch: float   # julian date
#     jdsatepochF: float  # julian date fraction
#     no_kozai: float     # kozai mean motion [rev/day], line 2, ch 53-63
#     ecco: float         # eccentricity, line 2, ch 27-33
#     sma: float          # semi-major axis [km]
#     inclo: float        # inclination [deg], line 2, ch 9-16
#     nodeo: float        # right ascension of ascending node [deg], line 2, ch 18-25
#     argpo: float        # argument of perigee [deg] line 2, ch 35-42
#     mo: float           # mean anomolay, line 2, ch 44-51
#     ndot: float         # first derivative of mean motion, line 1 [rad/s]
#     nddot: float        # second derivative of mean motion, line 1  [rad/s^2]
#     ndot_raw: float     # first derivative of mean motion divided by 2, line 1 [rev/day^2]
#     nddot_raw: float    # second derivative of mean motion divided by 6, line 1  [rev/day^3]
#     bstar: float        # B star drag term, line 1
#     elnum: int          # element number, line 1
#     revnum: int         # revolution number at epoch
#     classification: str = 'U'
#     ephtype: int = 0    # element set type

#     @classmethod
#     def from_tle(cls, tle1, tle2):
#         return tle_to_omm(tle1, tle2)

#     @property
#     def epoch(self):
#         """
#         Return julian date epoch
#         """
#         return self.jdsatepoch + self.jdsatepochF


# def tle_to_omm(tle1: str, tle2: str) -> OMM:
#     """
#     Convert TLE strings to OMM data
#     """
#     satnum = tle1[2:7]
#     classification = tle1[8]
#     epoch_year = int(tle1[18:20])
#     epoch_days = float(tle1[20:32])
#     jdsatepoch, jdsatepochF = epoch_to_jd(epoch_year, epoch_days)
#     # convert derivative of motion from rev/day to rad/s
#     ndot_raw = float(tle1[33:44])
#     nddot_raw = float(tle1[45:50]) * (10 ** float(tle1[50:52]))
#     ndot = ndot_raw * 2* (2*math.pi/(86400.0**2))
#     nddot = nddot_raw * 6 * (2*math.pi/(86400.0**3))
#     bstar = float(tle1[54:59]) * (10 ** float(tle1[59:61]))
#     ephtype = tle1[63]
#     elnum = int(tle1[65:69])

#     inclo = float(tle2[9:17])  # inclination
#     nodeo = float(tle2[18:26])  # right ascension of ascending node
#     ecco = float(tle2[27:34]) / 1e7  # eccentricity
#     argpo = float(tle2[35:43])
#     mo = float(tle2[43:52])    # mean anomaly
#     no_kozai = float(tle2[53:64])   # mean motion
#     sma = 6.6228 / (no_kozai**(2/3))  # ref: http://sat.belastro.net/satelliteorbitdetermination.com/orbit_elements_wiki.htm
#     revnum = int(tle2[64:69])

#     omm = OMM(
#         jdsatepoch=jdsatepoch,
#         jdsatepochF=jdsatepochF,
#         ndot_raw=ndot_raw,
#         nddot_raw=nddot_raw,
#         ndot=ndot,
#         nddot=nddot,
#         bstar=bstar,
#         inclo=inclo,
#         nodeo=nodeo,
#         ecco=ecco,
#         argpo=argpo,
#         mo=mo,
#         no_kozai=no_kozai,
#         revnum=revnum,
#         elnum=elnum,
#         classification=classification,
#         ephtype=ephtype
#     )

#     return omm


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