from __future__ import annotations
from math import pow, sin, cos, asin, sqrt, tan, pi, atan, degrees
from datetime import datetime, timedelta, timezone
from typing import NamedTuple, Union, Tuple

import numpy as np
from sgp4.earth_gravity import wgs84
from orbit_predictor.keplerian import rv2coe

from passpredict._time import jday2datetime_us, epoch_to_jd
from passpredict.constants import MU_E



def nu_to_anomaly(nu: float, e: float) -> float:
    """
    Find Anomaly. For elliptical orbits only
    Ref: Vallado, Alg 5, p.77
    """
    if e < 1:
        E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(0.5*nu))
    else:
        raise Exception('Orbit is not elliptical')
    return E


def anomaly_to_nu(E: float , e: float) -> float:
    """
    For elliptical orbits only
    Ref: Vallado, Alg 6, p.77
    """
    if 0 <= e < 1:
        nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(0.5*E))
    else:
        raise Exception('Orbit is not elliptical')
    return nu


def kepler_eqn_E(M: float, e: float) -> float:
    """
    Solve the Kepler equation for elliptical orbits

    Ref: Vallado, Alg 2, p.65
    """
    if (-pi < M < 0) or (M > pi):
        E = M - e
    else:
        E = M + e
    tol = 1e-8
    err = 1
    # Newton Raphson iteration
    while err > tol:
        E_next = E + (M - E + e*sin(E)) / (1 - e*cos(E))
        err = abs(E_next - E)
        E = E_next
    return E


def kepler_eqn(E: float, e: float) -> float:
    """ Get mean anomaly from eccentric anomaly """
    return E - e*sin(E)


def unkozai(no_kozai, ecco, inclo, whichconst):
    """
    Undo Kozai transformation
    Ref: orbit_predictor/utils.py, Vallado p.693, Sec 9.7.2
    """
    _, _, _, xke, j2, _, _, _ = whichconst
    ak = pow(xke / no_kozai, 2.0 / 3.0)
    d1 = 0.75 * j2 * (3.0 * cos(inclo)**2 - 1.0) / (1.0 - ecco**2)**(3/2)
    del_ = d1 / (ak * ak)
    adel = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
    return no_kozai / (1.0 + d1/(adel*adel))


class Orbit:
    def __init__(self, *a, **kw):
        self.satid = 0
        self.name = ""
        self._jdepoch = 0      # julian date
        self._jdepochF = 0     # julian date fraction
        self._no_kozai = 0     # kozai mean motion [rev/day], line 2, ch 53-63
        self._n = 0            # mean motion [rev/day]
        self._ecc = 0          # eccentricity, line 2, ch 27-33
        self._a = 0          # semi-major axis [km]
        self._inc = 0          # inclination [deg], line 2, ch 9-16
        self._raan = 0         # right ascension of ascending node [deg], line 2, ch 18-25
        self._argp = 0         # argument of perigee [deg] line 2, ch 35-42
        self._M = 0            # mean anomolay, line 2, ch 44-51
        self._nu = 0           # true anomaly
        self._ndot = 0         # first derivative of mean motion, line 1 [rad/s]
        self._nddot = 0        # second derivative of mean motion, line 1  [rad/s^2]
        self._ndot_raw = 0     # first derivative of mean motion divided by 2, line 1 [rev/day^2]
        self._nddot_raw = 0    # second derivative of mean motion divided by 6, line 1  [rev/day^3]
        self._bstar = 0        # B star drag term, line 1

    @property
    def epoch(self):
        jd = self._jdepoch + self._jdepochF
        return jday2datetime_us(jd)

    @property
    def no_kozai(self):
        """ Kozai mean motion [rev/day] """
        return self._no_kozai

    @property
    def n(self):
        """ Mean motion [rev/day] """
        return self._n

    @property
    def jdepoch(self):
        return self._jdepoch

    @property
    def jdepochF(self):
        return self._jdepochF

    @property
    def ecc(self):
        return self._ecc

    @property
    def bstar(self):
        return self._bstar

    @property
    def argp(self):
        return self._argp

    @property
    def a(self):
        return self._a

    @property
    def inc(self):
        return self._inc

    @property
    def nu(self):
        """ True anomaly [deg] """
        return self._nu

    @property
    def M(self):
        """ Mean anomaly [deg] """
        return self._M

    @property
    def raan(self):
        return self._raan

    @property
    def ndot(self):
        return self._ndot

    @property
    def nddot(self):
        return self._nddot

    @classmethod
    def from_tle(cls, tle: TLE):
        orbit = cls()
        orbit.satid = tle.satid
        orbit.name = tle.name
        orbit._jdepoch = tle.jdepoch
        orbit._bstar = tle.bstar
        orbit._ndot = tle.ndot * 2
        orbit._nddot = tle.nddot * 6
        orbit._inc = tle.inc
        orbit._ecc = tle.ecc
        orbit._raan = tle.raan
        orbit._M = tle.M
        orbit._argp = tle.argp

        # mean motion
        orbit._no_kozai = tle.no_kozai
        orbit._n = unkozai(tle.no_kozai, tle.ecc, tle.inc, wgs84)
        a = pow(MU_E / (tle.no_kozai * tle.no_kozai), 1/3)
        orbit._a = a   # semi major axis
        return orbit

    @classmethod
    def from_coe(cls, n, ecc, inc, argp, raan, M, jdepoch=0, ndot=0, nddot=0):
        """
        Create Orbit object from classical orbital elements
        """
        orbit = cls()
        orbit._inc = inc
        orbit._ecc = ecc
        orbit._ndot = ndot
        orbit._nddot = nddot
        jd, jdfr = divmod(jdepoch, 1)
        orbit._jdepoch = jd
        orbit._jdepochF = jdfr
        orbit._argp = argp
        orbit._raan = raan
        # mean motion
        orbit._n = n          # mean motion
        orbit._no_kozai = 0   # kozai mean motion
        # anomalies
        E = kepler_eqn_E(M, ecc)   # get eccentric anomaly from Kepler equation
        nu = anomaly_to_nu(E, ecc)  # then get the true anomaly
        orbit._M = M   # mean anomaly
        orbit._nu = nu  # true anomaly
        return orbit

    @classmethod
    def from_rv(cls, r, v, jdepoch=0, ndot=0, nddot=0):
        """
        Create Orbit object from ECI position and velocity vectors
        """
        orbit = cls()
        p, ecc, inc, raan, argp, nu = rv2coe(MU_E, r, v)
        orbit._inc = degrees(inc)
        orbit._ecc = ecc
        orbit._argp = degrees(argp)
        orbit._raan = degrees(raan)
        # semi parameter
        a = p/(1-ecc*ecc)
        orbit._a = a
        # mean motion
        n = sqrt(MU_E / (a*a*a))
        orbit._n = n          # mean motion
        # anomalies
        orbit._nu = nu  # true anomaly
        # derivative of mean motion
        orbit._ndot = ndot
        orbit._nddot = nddot
        # epoch
        jd, jdfr = divmod(jdepoch, 1)
        orbit._jdepoch = jd
        orbit._jdepochF = jdfr
        return orbit


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
        return epoch_from_tle(self.tle1)

    @property
    def tle1(self) -> str:
        return self.lines[0]

    @property
    def tle2(self) -> str:
        return self.lines[1]

    @property
    def intldesg(self):
        return self.tle1[9:17].rstrip()

    @property
    def jdepoch(self):
        epoch_string = self.tle1[18:32]
        year = int(epoch_string[0:2])
        days = float(epoch_string[2:])
        jd, jdfr = epoch_to_jd(year, days)
        jd = jd + jdfr
        return jd

    @property
    def ndot(self):
        # ndot = ndot_raw * 2* (2*math.pi/(86400.0**2))
        return float(self.tle1[33:44])

    @property
    def nddot(self):
        # nddot = nddot_raw * 6 * (2*math.pi/(86400.0**3))
        return float(self.tle1[45:50]) * pow(10,float(self.tle1[50:52]))

    @property
    def bstar(self):
        res = float(self.tle1[53] + '.' + self.tle1[54:59])
        return res * pow(10, float(self.tle1[59:61]))

    @property
    def ephtype(self):
        return float(self.tle1[62])

    @property
    def elnum(self):
        return float(self.tle1[65:69])

    @property
    def inc(self):
        """ Inclination [deg] """
        return float(self.tle2[9:17])  # inclination

    @property
    def raan(self):
        """ Right ascension of ascending node [deg] """
        return float(self.tle2[17:25])  # right ascension of ascending node

    @property
    def ecc(self):
        """ Eccentricity """
        return float('0.' + self.tle2[26:33].replace(' ','0'))  # eccentricity

    @property
    def argp(self):
        """ Argument of perigee [deg] """
        return float(self.tle2[34:42])

    @property
    def M(self):
        """ Mean anomaly [deg] """
        return float(self.tle2[43:52])    # mean anomaly

    @property
    def no_kozai(self):
        """ Kozai mean motion [rev/day] """
        return float(self.tle2[52:63])   # kozai mean motion

    # @cached_property
    # def sma(self):
    #     return 6.6228 / (self.no_kozai**(2/3))  # ref: http://sat.belastro.net/satelliteorbitdetermination.com/orbit_elements_wiki.htm

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