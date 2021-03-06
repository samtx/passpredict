import math
from math import pi

import numpy as np

from .nutation import fk5_nutation
from .precession import fk5_precession
from .core import rot3
from passpredict.constants import tau, ASEC2RAD, DEG2RAD
from passpredict.topocentric import site_ECEF, site_declination_and_K
from passpredict.timefn import jd2jc

# _arrays = np.load('nutation.npz')
# lunisolar_longitude_coefficients = _arrays['lunisolar_longitude_coefficients']
# lunisolar_obliquity_coefficients = _arrays['lunisolar_obliquity_coefficients']
# nals_t = _arrays['nals_t']
# napl_t = _arrays['napl_t']
# nutation_coefficients_longitude = _arrays['nutation_coefficients_longitude']
# nutation_coefficients_obliquity = _arrays['nutation_coefficients_obliquity']

# @profile



def site2eci(lat, lon, h, jdt):
    """Compute site ECI vector for jdt array

    TO DO: write this function

    """
    rsiteECI = 0
    return rsiteECI


def fk5(r, xp=0.0, yp=0.0):
    """IAU-76 / FK5 reductions for polar motion, nutation, precession

    Args
        r : position vector in ECI (ITRF) coordinates
        xp : polar motion along x axis in radians
        yp : polar motion along y axis in radians
    References:
        Vallado, Alg. 24, p.228
    """

    # polar motion using small angle approximations
    # Ref: Vallado, Eq 3-78
    rW = np.empty(r.shape)
    rW[0] = (1 + xp) * r[0]
    rW[1] = (1 - yp) * r[1]
    rW[2] = (1 - xp + yp) * r[2]

    # nutation, from IAU-1980 Theory of Nutation

    # We need to find the precession and nutation angles:
    # z, Theta, zeta, epsilon, deltaPsi1980, epsbar1980, thetaGAST1982


def eps1982(tt):
    """Compute the average eps 1982

    Args:
        tt : terrestial time

    References:
        Vallado, p. 225, Eq. 3-81
    """
    return 23.439291 * DEG2RAD + (-0.0130042 + (-1.64e-7 + 5.04e-7 * tt) * tt) * tt


def omega_moon(tt):
    """Nutation parameters for the moon

    Args:
        tt : terrestial time

    References:
        Vallado, p. 225, Eq. 3-82
    """
    r = 360 * 60
    omega_moon_deg = 125.04455501
    omega_moon_deg += (-5 * r - 134.1361851 + (0.0020756 + 2.139e-6 * tt) * tt) * tt
    return omega_moon_deg * DEG2RAD


# def fk5_nutation(tt):
#     """
#     Ref: Table D-6, p. 1043
#     """
#     api = an1 * M_moon + an2 * M_sun + an3 * Um_moon + an4 * D_sun + an5 * Omega_moon
#     psi_tmp = np.sum(A + B * tt) * np.sin(api)
#     eps_tmp = np.sum(C + D * tt) * np.cos(api)

#     return dPsi1980, dEps1980



    # # Convert PEF to GCRS with IAU2000a conventions

    # # We want to use IAU2000a nutation with IAU2006 precession

    # # we want to convert ITRF to GCRF

    # # rotate z by (t_gast*2pi/24) radians

    # # Find equation of equinoxes
    # dp, de = iau2000a(tt) # nutation angles
    # c_terms = equation_of_the_equinoxes_complimentary_terms(t.tt) / ASEC2RAD

    # d_psi = dp * 1e-7 + t.psi_correction
    # d_eps = de * 1e-7 + t.eps_correction

    # mean_ob = (mean_obliquity(t.tdb))/3600.0
    # true_ob = mean_ob + d_eps

    # eq_equinox = (d_psi * np.cos(mean_ob * DEG2RAD) + c_terms)/15.0

    # equation_of_equinoxes = earth_tilt[2] # [seconds of time]

    # t_gast = t_gmst + equation_of_equinoxes/3600


    # return rECEF


def iau2000a(jd_tt):
    """Compute Earth nutation based on the IAU 2000A nutation model.

    `jd_tt` - Terrestrial Time: Julian date float, or NumPy array of floats

    Returns a tuple ``(delta_psi, delta_epsilon)`` measured in tenths of
    a micro-arcsecond.  Each value is either a float, or a NumPy array
    with the same dimensions as the input argument.

    """
    # Interval between fundamental epoch J2000.0 and given date.

    t = (jd_tt - T0) / 36525.0

    # Compute fundamental arguments from Simon et al. (1994), in radians.

    a = fundamental_arguments(t)

    # ** Luni-solar nutation **
    # Summation of luni-solar nutation series (in reverse order).

    arg = nals_t.dot(a)
    np.fmod(arg, tau, out=arg)

    sarg = np.sin(arg)
    carg = np.cos(arg)

    stsc = np.array((sarg, t * sarg, carg)).T
    ctcs = np.array((carg, t * carg, sarg)).T

    dpsi = np.tensordot(stsc, lunisolar_longitude_coefficients)
    deps = np.tensordot(ctcs, lunisolar_obliquity_coefficients)

    # Compute and add in planetary components.

    if getattr(t, 'shape', ()) == ():
        a = t * anomaly_coefficient + anomaly_constant
    else:
        a = (np.outer(anomaly_coefficient, t).T + anomaly_constant).T
    a[-1] *= t

    np.fmod(a, tau, out=a)
    arg = napl_t.dot(a)
    np.fmod(arg, tau, out=arg)
    sc = np.array((sin(arg), cos(arg))).T

    dpsi += np.tensordot(sc, nutation_coefficients_longitude)
    deps += np.tensordot(sc, nutation_coefficients_obliquity)

    return dpsi, deps


def equation_of_the_equinoxes_complimentary_terms(jd_tt):
    """Compute the complementary terms of the equation of the equinoxes.

    `jd_tt` - Terrestrial Time: Julian date float, or NumPy array of floats

    """
    # Interval between fundamental epoch J2000.0 and current date.

    t = (jd_tt - T0) / 36525.0

    # Build array for intermediate results.

    shape = getattr(jd_tt, 'shape', ())
    fa = zeros((14,) if shape == () else (14, shape[0]))

    # Mean Anomaly of the Moon.

    fa[0] = ((485868.249036 +
              (715923.2178 +
              (    31.8792 +
              (     0.051635 +
              (    -0.00024470)
              * t) * t) * t) * t) * ASEC2RAD
              + (1325.0*t % 1.0) * tau)

    # Mean Anomaly of the Sun.

    fa[1] = ((1287104.793048 +
              (1292581.0481 +
              (     -0.5532 +
              (     +0.000136 +
              (     -0.00001149)
              * t) * t) * t) * t) * ASEC2RAD
              + (99.0*t % 1.0) * tau)

    # Mean Longitude of the Moon minus Mean Longitude of the Ascending
    # Node of the Moon.

    fa[2] = (( 335779.526232 +
              ( 295262.8478 +
              (    -12.7512 +
              (     -0.001037 +
              (      0.00000417)
              * t) * t) * t) * t) * ASEC2RAD
              + (1342.0*t % 1.0) * tau)

    # Mean Elongation of the Moon from the Sun.

    fa[3] = ((1072260.703692 +
              (1105601.2090 +
              (     -6.3706 +
              (      0.006593 +
              (     -0.00003169)
              * t) * t) * t) * t) * ASEC2RAD
              + (1236.0*t % 1.0) * tau)

    # Mean Longitude of the Ascending Node of the Moon.

    fa[4] = (( 450160.398036 +
              (-482890.5431 +
              (      7.4722 +
              (      0.007702 +
              (     -0.00005939)
              * t) * t) * t) * t) * ASEC2RAD
              + (-5.0*t % 1.0) * tau)

    fa[ 5] = (4.402608842 + 2608.7903141574 * t)
    fa[ 6] = (3.176146697 + 1021.3285546211 * t)
    fa[ 7] = (1.753470314 +  628.3075849991 * t)
    fa[ 8] = (6.203480913 +  334.0612426700 * t)
    fa[ 9] = (0.599546497 +   52.9690962641 * t)
    fa[10] = (0.874016757 +   21.3299104960 * t)
    fa[11] = (5.481293872 +    7.4781598567 * t)
    fa[12] = (5.311886287 +    3.8133035638 * t)
    fa[13] = (0.024381750 +    0.00000538691 * t) * t

    fa %= tau

    # Evaluate the complementary terms.

    a = ke0_t.dot(fa)
    s0 = se0_t_0.dot(sin(a)) + se0_t_1.dot(cos(a))

    a = ke1.dot(fa)
    s1 = se1_0 * sin(a) + se1_1 * cos(a)

    c_terms = s0 + s1 * t
    c_terms *= ASEC2RAD
    return c_terms

fa0, fa1, fa2, fa3, fa4 = np.array((

    # Mean Anomaly of the Moon.
    (485868.249036, 1717915923.2178, 31.8792, 0.051635, - .00024470),

    # Mean Anomaly of the Sun.
    (1287104.79305,  129596581.0481, - 0.5532, 0.000136, - 0.00001149),

    # Mean Longitude of the Moon minus Mean Longitude of the Ascending
    # Node of the Moon.
    (335779.526232, 1739527262.8478, - 12.7512, -  0.001037, 0.00000417),

    # Mean Elongation of the Moon from the Sun.
    (1072260.70369, 1602961601.2090, - 6.3706, 0.006593, - 0.00003169),

    # Mean Longitude of the Ascending Node of the Moon.
    (450160.398036, - 6962890.5431, 7.4722, 0.007702, - 0.00005939),

    )).T[:,:,None]

def fundamental_arguments(t):

    """Compute the fundamental arguments (mean elements) of Sun and Moon.

    `t` - TDB time in Julian centuries since J2000.0, as float or NumPy array

    Outputs fundamental arguments, in radians:
          a[0] = l (mean anomaly of the Moon)
          a[1] = l' (mean anomaly of the Sun)
          a[2] = F (mean argument of the latitude of the Moon)
          a[3] = D (mean elongation of the Moon from the Sun)
          a[4] = Omega (mean longitude of the Moon's ascending node);
                 from Simon section 3.4(b.3),
                 precession = 5028.8200 arcsec/cy)

    """
    a = fa4 * t
    a += fa3
    a *= t
    a += fa2
    a *= t
    a += fa1
    a *= t
    a += fa0
    fmod(a, ASEC360, out=a)
    a *= ASEC2RAD
    if getattr(t, 'shape', ()):
        return a
    return a[:,0]






# @profile



##################
# From sgp4lib.py in skyfield
###################


def ITRF_position_velocity_error(t):
    """Return the ITRF position, velocity, and error at time `t`.

    The position is an x,y,z vector measured in au, the velocity is
    an x,y,z vector measured in au/day, and the error is a vector of
    possible error messages for the time or vector of times `t`.

    """
    rTEME, vTEME, error = sgp4(tle, t)
    rTEME /= AU_KM
    vTEME /= AU_KM
    vTEME *= DAY_S
    rITRF, vITRF = TEME_to_ITRF(t.ut1, rTEME, vTEME)
    return rITRF, vITRF, error


def _at(self, t):
    """Compute this satellite's GCRS position and velocity at time `t`."""
    rITRF, vITRF, error = self.ITRF_position_velocity_error(t)
    rGCRS, vGCRS = ITRF_to_GCRS2(t, rITRF, vITRF)
    return rGCRS, vGCRS, rGCRS, error


def TEME_to_ITRF(jd_ut1, rTEME, vTEME, xp=0.0, yp=0.0):
    """Convert TEME position and velocity into standard ITRS coordinates.

    This converts a position and velocity vector in the idiosyncratic
    True Equator Mean Equinox (TEME) frame of reference used by the SGP4
    theory into vectors into the more standard ITRS frame of reference.
    The velocity should be provided in units per day, not per second.

    From AIAA 2006-6753 Appendix C.
    """
    theta, theta_dot = theta_GMST1982(jd_ut1)
    zero = theta_dot * 0.0
    angular_velocity = np.array([zero, zero, -theta_dot])
    R = rot3(-theta).T
    if len(rTEME.shape) == 1:
        rPEF = np.dot(R, rTEME)
        vPEF = np.dot(R, vTEME) + np.cross(angular_velocity, rPEF)
    else:
        rPEF = np.einsum("ij...,j...->i...", R, rTEME)
        vPEF = (
            np.einsum("ij...,j...->i...", R, vTEME)
            + np.cross(angular_velocity, rPEF, 0, 0).T
        )
    if xp == 0.0 and yp == 0.0:
        rITRF = rPEF
        vITRF = vPEF
    else:
        W = (rot1(yp)).dot(rot2(xp))
        W = W.T
        rITRF = (W).dot(rPEF)
        vITRF = (W).dot(vPEF)
    return rITRF, vITRF


def ITRF_to_GCRS2(t, rITRF, vITRF):
    """
    From skyfield
    url: https://github.com/skyfielders/python-skyfield/blob/5a9c1a0f5254093d95ea66cb12a9498f24440326/skyfield/positionlib.py#L691

    """

    # Todo: wobble

    spin = rot_z(t.gast * tau / 24.0)

    position = np.einsum('ij...,j...->i...', spin, array(rITRF))
    position = np.einsum('ij...,j...->i...', t.MT, position)

    velocity = np.einsum('ij...,j...->i...', spin, array(vITRF))
    velocity = np.einsum('ij...,j...->i...', t.MT, velocity)
    velocity[0] += DAY_S * ANGVEL * - position[1]
    velocity[1] += DAY_S * ANGVEL * position[0]

    return position, velocity


## from timelib.py in skyfield

def gmst(self):
    return sidereal_time(self)

def gast(self):
    return self.gmst + self._earth_tilt[2] / 3600.0

def _earth_tilt(self):
    return earth_tilt(self)