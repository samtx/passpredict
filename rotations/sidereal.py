import numpy as np

from passpredict.constants import ASEC2RAD, DEG2RAD, RAD2DEG, tau


def equinox1982(dPsi1980, eps1980, omega_moon):
    """Return the equation of equinoxes for sidereal time
    for IAU-76/FK5 reductions

    Args:
        t : time

    References:
        Vallado, p.224
    """
    eq = equinox1982_geometric_terms(dPsi1980, eps1980)
    eq += 0.00264 * ASEC2RAD * np.sin(omega_moon)
    eq += 0.000063 * np.sin(2 * omega_moon)
    return eq


def equinox1982_geometric_terms(dPsi1980, eps1980):
    """Return the geometric terms of equation of equinoxes for sidereal time
    for IAU-76/FK5 reductions

    Args:
        dPsi1980 : float  (radians)
        eps1920 : float   (radians)

    References:
        Vallado, p.224
    """
    return dPsi1980 * np.cos(eps1980)


def theta_GAST1982(Eq, gmst):
    """Return the Greenwich apparent sidereal time

    Sidereal time for IAU-76/FK5 reductions

    Args:
        Eq : equinox1982
        gmst : Greenwich mean sidreal time

    References:
        Vallado, p.224
    """
    return Eq + gmst


def theta_GMST1982(jd_ut1):
    """Return the angle of Greenwich Mean Standard Time 1982 given the JD.

    This angle defines the difference between the idiosyncratic True
    Equator Mean Equinox (TEME) frame of reference used by SGP4 and the
    more standard Pseudo Earth Fixed (PEF) frame of reference.

    Args:
        jd_ut1 (float): Julian date since Jan 1, 2000

    Returns:
        theta (float): GMST in radians
        thetadt (float): time derivative of GMST in rad/s

    References:
        Vallado, et al. "Revisiting Spacetrack Report #3", AIAA, 2006-6753, Appendix C.
        Rhodes, Skyfield library, github.com/skyfielders/python-skyfield --> sgp4lib.py
    """
    tau = 2 * np.pi
    _second = 1.0 / (24.0 * 60.0 * 60.0)
    T0 = 2451545.0  # JD for Jan 1, 2000
    t = (jd_ut1 - T0) / 36525.0
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t
    dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t
    theta = (jd_ut1 % 1.0 + g * _second % 1.0) * tau
    theta_dot = (1.0 + dg * _second / 36525.0) * tau
    return theta, theta_dot