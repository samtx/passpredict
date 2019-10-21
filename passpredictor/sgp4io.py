"""Read the TLE earth satellite file format.

This is a minimally-edited copy of "sgp4io.cpp".

"""
import re
from datetime import datetime
from math import pi, pow
from passpredictor.sgp4 import sgp4init
from collections import namedtuple
from passpredictor.sgp4 import getgravconst

INT_RE = re.compile(r'[+-]?\d*')
FLOAT_RE = re.compile(r'[+-]?\d*(\.\d*)?')

LINE1 = '1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN'
LINE2 = '2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN'

error_message = """TLE format error

The Two-Line Element (TLE) format was designed for punch cards, and so
is very strict about the position of every period, space, and digit.
Your line does not quite match.  Here is the official format for line {0}
with an N where each digit should go, followed by the line you provided:

{1}
{2}"""

"""
/*     ----------------------------------------------------------------
*
*                               sgp4io.cpp
*
*    this file contains a function to read two line element sets. while
*    not formerly part of the sgp4 mathematical theory, it is
*    required for practical implemenation.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              27 Aug 10  david vallado
*                           fix input format and delete unused variables in twoline2rv
*    changes :
*               3 sep 08  david vallado
*                           add operationmode for afspc (a) or improved (i)
*               9 may 07  david vallado
*                           fix year correction to 57
*              27 mar 07  david vallado
*                           misc fixes to manual inputs
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */
"""

"""Three earth-gravity models for use with SGP4."""
EarthGravity = namedtuple(
    'EarthGravity',
    'tumin mu radiusearthkm xke j2 j3 j4 j3oj2',
    )

wgs72old = EarthGravity(*getgravconst('wgs72old'))
wgs72 = EarthGravity(*getgravconst('wgs72'))
wgs84 = EarthGravity(*getgravconst('wgs84'))


"""Utility routines from "sgp4ext.cpp"."""

from math import (acos, asinh, atan2, copysign, cos, fabs, fmod,
                  pi, sin, sinh, sqrt, tan)

undefined = None

"""
/* -----------------------------------------------------------------------------
*
*                           function mag
*
*  this procedure finds the magnitude of a vector.  the tolerance is set to
*    0.000001, thus the 1.0e-12 for the squared test of underflows.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec         - vector
*
*  outputs       :
*    vec         - answer stored in fourth component
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */
"""

def mag(x):
     return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

"""
/* -----------------------------------------------------------------------------
*
*                           procedure cross
*
*  this procedure crosses two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    outvec      - vector result of a x b
*
*  locals        :
*    none.
*
*  coupling      :
*    mag           magnitude of a vector
 ---------------------------------------------------------------------------- */
"""

def cross(vec1, vec2, outvec):
     outvec[0]= vec1[1]*vec2[2] - vec1[2]*vec2[1];
     outvec[1]= vec1[2]*vec2[0] - vec1[0]*vec2[2];
     outvec[2]= vec1[0]*vec2[1] - vec1[1]*vec2[0];

"""
/* -----------------------------------------------------------------------------
*
*                           function dot
*
*  this function finds the dot product of two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    dot         - result
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
* --------------------------------------------------------------------------- */
"""

def dot(x, y):
     return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);

"""
/* -----------------------------------------------------------------------------
*
*                           procedure angle
*
*  this procedure calculates the angle between two vectors.  the output is
*    set to 999999.1 to indicate an undefined value.  be sure to check for
*    this at the output phase.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    theta       - angle between the two vectors  -pi to pi
*
*  locals        :
*    temp        - temporary real variable
*
*  coupling      :
*    dot           dot product of two vectors
* --------------------------------------------------------------------------- */
"""

def angle(vec1, vec2):

     small     = 0.00000001;
     undefined = 999999.1;

     magv1 = mag(vec1);
     magv2 = mag(vec2);

     if magv1*magv2 > small*small:

         temp= dot(vec1,vec2) / (magv1*magv2);
         if fabs(temp) > 1.0:
             temp = copysign(1.0, temp)
         return acos( temp );

     else:
         return undefined;

"""
/* -----------------------------------------------------------------------------
*
*                           function newtonnu
*
*  this function solves keplers equation when the true anomaly is known.
*    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
*    the parabolic limit at 168° is arbitrary. the hyperbolic anomaly is also
*    limited. the hyperbolic sine is used because it's not double valued.
*
*  author        : david vallado                  719-573-2600   27 may 2002
*
*  revisions
*    vallado     - fix small                                     24 sep 2002
*
*  inputs          description                    range / units
*    ecc         - eccentricity                   0.0  to
*    nu          - true anomaly                   -2pi to 2pi rad
*
*  outputs       :
*    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 °
*    m           - mean anomaly                   0.0  to 2pi rad       151.7425 °
*
*  locals        :
*    e1          - eccentric anomaly, next value  rad
*    sine        - sine of e
*    cose        - cosine of e
*    ktr         - index
*
*  coupling      :
*    asinh       - arc hyperbolic sine
*
*  references    :
*    vallado       2007, 85, alg 5
* --------------------------------------------------------------------------- */
"""

def newtonnu(ecc, nu):

     #  ---------------------  implementation   ---------------------
     e0= 999999.9;
     m = 999999.9;
     small = 0.00000001;

     #  --------------------------- circular ------------------------
     if fabs(ecc) < small:

         m = nu;
         e0= nu;

     else:
         #  ---------------------- elliptical -----------------------
         if ecc < 1.0-small:

             sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
             cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
             e0  = atan2( sine,cose );
             m   = e0 - ecc*sin(e0);

         else:
             #  -------------------- hyperbolic  --------------------
             if ecc > 1.0 + small:

                 if ecc > 1.0 and fabs(nu)+0.00001 < pi-acos(1.0 /ecc):

                     sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                     e0  = asinh( sine );
                     m   = ecc*sinh(e0) - e0;

             else:
                 #  ----------------- parabolic ---------------------
                 if fabs(nu) < 168.0*pi/180.0:

                     e0= tan( nu*0.5  );
                     m = e0 + (e0*e0*e0)/3.0;

     if ecc < 1.0:

         m = fmod( m,2.0 *pi );
         if m < 0.0:
             m = m + 2.0 *pi;
         e0 = fmod( e0,2.0 *pi );

     return e0, m


"""
/* -----------------------------------------------------------------------------
*
*                           function rv2coe
*
*  this function finds the classical orbital elements given the geocentric
*    equatorial position and velocity vectors.
*
*  author        : david vallado                  719-573-2600   21 jun 2002
*
*  revisions
*    vallado     - fix special cases                              5 sep 2002
*    vallado     - delete extra check in inclination code        16 oct 2002
*    vallado     - add constant file use                         29 jun 2003
*    vallado     - add mu                                         2 apr 2007
*
*  inputs          description                    range / units
*    r           - ijk position vector            km
*    v           - ijk velocity vector            km / s
*    mu          - gravitational parameter        km3 / s2
*
*  outputs       :
*    p           - semilatus rectum               km
*    a           - semimajor axis                 km
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
*
*  locals        :
*    hbar        - angular momentum h vector      km2 / s
*    ebar        - eccentricity     e vector
*    nbar        - line of nodes    n vector
*    c1          - v**2 - u/r
*    rdotv       - r dot v
*    hk          - hk unit vector
*    sme         - specfic mechanical energy      km2 / s2
*    i           - index
*    e           - eccentric, parabolic,
*                  hyperbolic anomaly             rad
*    temp        - temporary variable
*    typeorbit   - type of orbit                  ee, ei, ce, ci
*
*  coupling      :
*    mag         - magnitude of a vector
*    cross       - cross product of two vectors
*    angle       - find the angle between two vectors
*    newtonnu    - find the mean anomaly
*
*  references    :
*    vallado       2007, 126, alg 9, ex 2-5
* --------------------------------------------------------------------------- */
"""

def rv2coe(r, v, mu):

     hbar = [None, None, None]
     nbar = [None, None, None]
     ebar = [None, None, None]
     typeorbit = [None, None, None];

     twopi  = 2.0 * pi;
     halfpi = 0.5 * pi;
     small  = 0.00000001;
     undefined = 999999.1;
     infinite  = 999999.9;

     #  -------------------------  implementation   -----------------
     magr = mag( r );
     magv = mag( v );

     #  ------------------  find h n and e vectors   ----------------
     cross( r,v, hbar );
     magh = mag( hbar );
     if magh > small:

         nbar[0]= -hbar[1];
         nbar[1]=  hbar[0];
         nbar[2]=   0.0;
         magn = mag( nbar );
         c1 = magv*magv - mu /magr;
         rdotv = dot( r,v );
         for i in range(0, 3):
             ebar[i]= (c1*r[i] - rdotv*v[i])/mu;
         ecc = mag( ebar );

         #  ------------  find a e and semi-latus rectum   ----------
         sme= ( magv*magv*0.5  ) - ( mu /magr );
         if fabs( sme ) > small:
             a= -mu  / (2.0 *sme);
         else:
             a= infinite;
         p = magh*magh/mu;

         #  -----------------  find inclination   -------------------
         hk= hbar[2]/magh;
         incl= acos( hk );

         #  --------  determine type of orbit for later use  --------
         #  ------ elliptical, parabolic, hyperbolic inclined -------
         typeorbit = 'ei'
         if ecc < small:

             #  ----------------  circular equatorial ---------------
             if  incl < small or fabs(incl-pi) < small:
                 typeorbit = 'ce'
             else:
                 #  --------------  circular inclined ---------------
                 typeorbit = 'ci'

         else:

             #  - elliptical, parabolic, hyperbolic equatorial --
             if incl < small or fabs(incl-pi) < small:
                 typeorbit = 'ee'

         #  ----------  find longitude of ascending node ------------
         if magn > small:

             temp= nbar[0] / magn;
             if fabs(temp) > 1.0:
                 temp = copysign(1.0, temp)
             omega= acos( temp );
             if nbar[1] < 0.0:
                 omega= twopi - omega;

         else:
             omega= undefined;

         #  ---------------- find argument of perigee ---------------
         if typeorbit == 'ei':

             argp = angle( nbar,ebar);
             if ebar[2] < 0.0:
                 argp= twopi - argp;

         else:
             argp= undefined;

         #  ------------  find true anomaly at epoch    -------------
         if typeorbit[0] == 'e':

             nu =  angle( ebar,r);
             if rdotv < 0.0:
                 nu= twopi - nu;

         else:
             nu= undefined;

         #  ----  find argument of latitude - circular inclined -----
         if typeorbit == 'ci':

             arglat = angle( nbar,r );
             if r[2] < 0.0:
                 arglat= twopi - arglat;
             m = arglat;

         else:
             arglat= undefined;

         #  -- find longitude of perigee - elliptical equatorial ----
         if ecc > small and typeorbit == 'ee':

             temp= ebar[0]/ecc;
             if fabs(temp) > 1.0:
                 temp = copysign(1.0, temp)
             lonper= acos( temp );
             if ebar[1] < 0.0:
                 lonper= twopi - lonper;
             if incl > halfpi:
                 lonper= twopi - lonper;

         else:
             lonper= undefined;

         #  -------- find true longitude - circular equatorial ------
         if magr > small and typeorbit == 'ce':

             temp= r[0]/magr;
             if fabs(temp) > 1.0:
                 temp = copysign(1.0, temp)
             truelon= acos( temp );
             if r[1] < 0.0:
                 truelon= twopi - truelon;
             if incl > halfpi:
                 truelon= twopi - truelon;
             m = truelon;

         else:
             truelon= undefined;

         #  ------------ find mean anomaly for all orbits -----------
         if typeorbit[0] == 'e':
             e, m = newtonnu(ecc, nu);

     else:
        p    = undefined;
        a    = undefined;
        ecc  = undefined;
        incl = undefined;
        omega= undefined;
        argp = undefined;
        nu   = undefined;
        m    = undefined;
        arglat = undefined;
        truelon= undefined;
        lonper = undefined;

     return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper

"""
/* -----------------------------------------------------------------------------
*
*                           procedure jday
*
*  this procedure finds the julian date given the year, month, day, and time.
*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
*
*  algorithm     : calculate the answer in one step for efficiency
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - universal time hour            0 .. 23
*    min         - universal time min             0 .. 59
*    sec         - universal time sec             0.0 .. 59.999
*
*  outputs       :
*    jd          - julian date                    days from 4713 bc
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 189, alg 14, ex 3-14
*
* --------------------------------------------------------------------------- */
"""

def jday(year, mon, day, hr, minute, sec):

  return (367.0 * year -
          7.0 * (year + ((mon + 9.0) // 12.0)) * 0.25 // 1.0 +
          275.0 * mon // 9.0 +
          day + 1721013.5 +
          ((sec / 60.0 + minute) / 60.0 + hr) / 24.0  #  ut in days
          #  - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
          )

"""
/* -----------------------------------------------------------------------------
*
*                           procedure days2mdhms
*
*  this procedure converts the day of the year, days, to the equivalent month
*    day, hour, minute and second.
*
*  algorithm     : set up array for the number of days per month
*                  find leap year - use 1900 because 2000 is a leap year
*                  loop through a temp value while the value is < the days
*                  perform int conversions to the correct day and month
*                  convert remainder into h m s using type conversions
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    days        - julian day of the year         0.0  .. 366.0
*
*  outputs       :
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - hour                           0 .. 23
*    min         - minute                         0 .. 59
*    sec         - second                         0.0 .. 59.999
*
*  locals        :
*    dayofyr     - day of year
*    temp        - temporary extended values
*    inttemp     - temporary int value
*    i           - index
*    lmonth[12]  - int array containing the number of days per month
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */
"""

def days2mdhms(year, days):

     lmonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

     dayofyr = int(days // 1.0);
     #  ----------------- find month and day of month ----------------
     if (year % 4) == 0:
       lmonth = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

     i = 1;
     inttemp = 0;
     while dayofyr > inttemp + lmonth[i-1] and i < 12:

       inttemp = inttemp + lmonth[i-1];
       i += 1;

     mon = i;
     day = dayofyr - inttemp;

     #  ----------------- find hours minutes and seconds -------------
     temp = (days - dayofyr) * 24.0;
     hr   = int(temp // 1.0);
     temp = (temp - hr) * 60.0;
     minute  = int(temp // 1.0);
     sec  = (temp - minute) * 60.0;

     return mon, day, hr, minute, sec

"""
/* -----------------------------------------------------------------------------
*
*                           procedure invjday
*
*  this procedure finds the year, month, day, hour, minute and second
*  given the julian date. tu can be ut1, tdt, tdb, etc.
*
*  algorithm     : set up starting values
*                  find leap year - use 1900 because 2000 is a leap year
*                  find the elapsed days through the year in a loop
*                  call routine to find each individual value
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    jd          - julian date                    days from 4713 bc
*
*  outputs       :
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - hour                           0 .. 23
*    min         - minute                         0 .. 59
*    sec         - second                         0.0 .. 59.999
*
*  locals        :
*    days        - day of year plus fractional
*                  portion of a day               days
*    tu          - julian centuries from 0 h
*                  jan 0, 1900
*    temp        - temporary double values
*    leapyrs     - number of leap years from 1900
*
*  coupling      :
*    days2mdhms  - finds month, day, hour, minute and second given days and year
*
*  references    :
*    vallado       2007, 208, alg 22, ex 3-13
* --------------------------------------------------------------------------- */
"""

def invjday(jd):

     #  --------------- find year and days of the year ---------------
     temp    = jd - 2415019.5;
     tu      = temp / 365.25;
     year    = 1900 + int(tu // 1.0);
     leapyrs = int(((year - 1901) * 0.25) // 1.0);

     #  optional nudge by 8.64x10-7 sec to get even outputs
     days    = temp - ((year - 1900) * 365.0 + leapyrs) + 0.00000000001;

     #  ------------ check for case of beginning of a year -----------
     if (days < 1.0):
         year    = year - 1;
         leapyrs = int(((year - 1901) * 0.25) // 1.0);
         days    = temp - ((year - 1900) * 365.0 + leapyrs);

     #  ----------------- find remaing data  -------------------------
     mon, day, hr, minute, sec = days2mdhms(year, days);
     sec = sec - 0.00000086400;
     return year, mon, day, hr, minute, sec


class Satellite(object):
    """An earth-orbiting satellite as represented by the SGP4 model.

    Most of this class's hundred-plus attributes are intermediate values
    of interest only to the propagation algorithm itself.  Here are the
    attributes set by ``sgp4.io.twoline2rv()`` in which users are likely
    to be interested:

    ``satnum``
        Unique satellite number given in the TLE file.
    ``epochyr``
        Full four-digit year of this element set's epoch moment.
    ``epochdays``
        Fractional days into the year of the epoch moment.
    ``jdsatepoch``
        Julian date of the epoch (computed from ``epochyr`` and ``epochdays``).
    ``ndot``
        First time derivative of the mean motion (ignored by SGP4).
    ``nddot``
        Second time derivative of the mean motion (ignored by SGP4).
    ``bstar``
        Ballistic drag coefficient B* in inverse earth radii.
    ``inclo``
        Inclination in radians.
    ``nodeo``
        Right ascension of ascending node in radians.
    ``ecco``
        Eccentricity.
    ``argpo``
        Argument of perigee in radians.
    ``mo``
        Mean anomaly in radians.
    ``no``
        Mean motion in radians per minute.

    """
    pass

"""
/* -----------------------------------------------------------------------------
*
*                           function twoline2rv
*
*  this function converts the two line element set character string data to
*    variables and initializes the sgp4 variables. several intermediate varaibles
*    and quantities are determined. note that the result is a structure so multiple
*    satellites can be processed simultaneously without having to reinitialize. the
*    verification mode is an important option that permits quick checks of any
*    changes to the underlying technical theory. this option works using a
*    modified tle file in which the start, stop, and delta time values are
*    included at the end of the second line of data. this only works with the
*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
*    from epoch and is useful when performing entire catalog runs.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs        :
*    longstr1    - first line of the tle
*    longstr2    - second line of the tle
*    typerun     - type of run                    verification 'v', catalog 'c',
*                                                 manual 'm'
*    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
*    opsmode     - mode of operation afspc or improved 'a', 'i'
*    whichconst  - which set of constants to use  72, 84
*
*  outputs       :
*    satrec      - structure containing all the sgp4 satellite information
*
*  coupling      :
*    getgravconst-
*    days2mdhms  - conversion of days to month, day, hour, minute, second
*    jday        - convert day month year hour minute second into julian date
*    sgp4init    - initialize the sgp4 variables
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */
"""

def twoline2rv(longstr1, longstr2, whichconst, opsmode='i'):
    """Return a Satellite imported from two lines of TLE data.

    Provide the two TLE lines as strings `longstr1` and `longstr2`,
    and select which standard set of gravitational constants you want
    by providing `gravity_constants`:

    `sgp4.earth_gravity.wgs72` - Standard WGS 72 model
    `sgp4.earth_gravity.wgs84` - More recent WGS 84 model
    `sgp4.earth_gravity.wgs72old` - Legacy support for old SGP4 behavior

    Normally, computations are made using various recent improvements
    to the algorithm.  If you want to turn some of these off and go
    back into "opsmode" mode, then set `opsmode` to `a`.

    """

    deg2rad  =   pi / 180.0;         #    0.0174532925199433
    xpdotp   =  1440.0 / (2.0 *pi);  #  229.1831180523293

    tumin = whichconst.tumin

    satrec = Satellite()
    satrec.error = 0;
    satrec.whichconst = whichconst  # Python extension: remembers its consts

    line = longstr1.rstrip()
    # try/except is not well supported by Numba
    if (len(line) >= 64 and
        line.startswith('1 ') and
        line[8] == ' ' and
        line[23] == '.' and
        line[32] == ' ' and
        line[34] == '.' and
        line[43] == ' ' and
        line[52] == ' ' and
        line[61] == ' ' and
        line[63] == ' '):

        _saved_satnum = satrec.satnum = int(line[2:7])
        satrec.line1 = line
        satrec.classification = line[7] or 'U'
        satrec.intldesg = line[9:17]
        two_digit_year = int(line[18:20])
        satrec.epochdays = float(line[20:32])
        satrec.ndot = float(line[33:43])
        satrec.nddot = float(line[44] + '.' + line[45:50])
        nexp = int(line[50:52])
        satrec.bstar = float(line[53] + '.' + line[54:59])
        ibexp = int(line[59:61])
        satrec.ephtype = line[62]
        satrec.elnum = int(line[64:68])
    else:
        raise ValueError(error_message.format(1, LINE1, line))

    line = longstr2.rstrip()
    if (len(line) >= 69 and
        line.startswith('2 ') and
        line[7] == ' ' and
        line[11] == '.' and
        line[16] == ' ' and
        line[20] == '.' and
        line[25] == ' ' and
        line[33] == ' ' and
        line[37] == '.' and
        line[42] == ' ' and
        line[46] == '.' and
        line[51] == ' '):

        satrec.satnum = int(line[2:7])
        if _saved_satnum != satrec.satnum:
            raise ValueError('Object numbers in lines 1 and 2 do not match')

        satrec.line2 = line
        satrec.inclo = float(line[8:16])
        satrec.nodeo = float(line[17:25])
        satrec.ecco = float('0.' + line[26:33].replace(' ', '0'))
        satrec.argpo = float(line[34:42])
        satrec.mo = float(line[43:51])
        satrec.no_kozai = float(line[52:63])
        satrec.revnum = line[63:68]
    #except (AssertionError, IndexError, ValueError):
    else:
        raise ValueError(error_message.format(2, LINE2, line))

    #  ---- find no, ndot, nddot ----
    satrec.no_kozai = satrec.no_kozai / xpdotp; #   rad/min
    satrec.nddot= satrec.nddot * pow(10.0, nexp);
    satrec.bstar= satrec.bstar * pow(10.0, ibexp);

    #  ---- convert to sgp4 units ----
    satrec.ndot = satrec.ndot  / (xpdotp*1440.0);  #   ? * minperday
    satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);

    #  ---- find standard orbital elements ----
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo  * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     * deg2rad;


    """
    // ----------------------------------------------------------------
    // find sgp4epoch time of element set
    // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    // and minutes from the epoch (time)
    // ----------------------------------------------------------------

    // ---------------- temp fix for years from 1957-2056 -------------------
    // --------- correct fix will occur when year is 4-digit in tle ---------
    """
    if two_digit_year < 57:
        year = two_digit_year + 2000;
    else:
        year = two_digit_year + 1900;

    mon,day,hr,minute,sec = days2mdhms(year, satrec.epochdays);
    sec_whole, sec_fraction = divmod(sec, 1.0)

    satrec.epochyr = year
    satrec.jdsatepoch = jday(year,mon,day,hr,minute,sec);
    satrec.epoch = datetime(year, mon, day, hr, minute, int(sec_whole),
                            int(sec_fraction * 1000000.0 // 1.0))

    #  ---------------- initialize the orbit at sgp4epoch -------------------
    sgp4init(whichconst, opsmode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar,
             satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo,
             satrec.no_kozai, satrec.nodeo, satrec)

    return satrec

def verify_checksum(*lines):
    """Verify the checksum of one or more TLE lines.

    Raises `ValueError` if any of the lines fails its checksum, and
    includes the failing line in the error message.

    """
    for line in lines:
        checksum = line[68:69]
        if not checksum.isdigit():
            continue
        checksum = int(checksum)
        computed = compute_checksum(line)
        if checksum != computed:
            complaint = ('TLE line gives its checksum as {}'
                         ' but in fact tallies to {}:\n{}')
            raise ValueError(complaint.format(checksum, computed, line))

def fix_checksum(line):
    """Return a new copy of the TLE `line`, with the correct checksum appended.

    This discards any existing checksum at the end of the line, if a
    checksum is already present.

    """
    return line[:68].ljust(68) + str(compute_checksum(line))

def compute_checksum(line):
    """Compute the TLE checksum for the given line."""
    return sum((int(c) if c.isdigit() else c == '-') for c in line[0:68]) % 10
