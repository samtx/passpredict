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


def get_overpasses(el, azm, rng, dt_ary, rSEZ, jdt=None, rsiteECI=None, rsatECI=None, loc=None, sat=None):
    el0 = el[:-1]
    el1 = el[1:]
    el_change_sign = (el0*el1 < 0)
    # Find the start of an overpass
    start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]
    # Find the end of an overpass
    end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]
    # Iterate over start/end indecies and gather inbetween indecies
    overpasses = np.empty(start_idx.size, dtype=object)
    for j in range(start_idx.size):
        # Store indecies of overpasses in a list
        idx0 = start_idx[j]
        idxf = end_idx[j]
        overpass_idx = np.arange(idx0, idxf+1, dtype=int)
        idxmax = np.argmax(el[overpass_idx])
        start_pt = Point(dt_ary[idx0], az[idx0], el[idx0], rng[idx0])
        max_pt = Point(dt_ary[idxmax], az[idxmax], el[idxmax], rng[idxmax])
        end_pt = Point(dt_ary[idxf], az[idxf], el[idxf], rng[idxf])
        sat_vis = satellite_visible(rsatECI, rsiteECI, rSEZ, jdt)
        overpass = Overpass(
            loc,
            sat,
            start_pt,
            max_pt,
            end_pt,
            dt_ary[overpass_idx],
            rSEZ[:,overpass_idx]
        )
        overpasses[j] = overpass
    return overpasses


def predict_passes(lat, lon, h, rsatECEF, rsatECI, jdt, rsun, loc=None, sat=None):
    rSEZ = site_sat_rotations(lat, lon, h, rsatECEF)
    rsiteECI = rotations.site2eci(lat, lon, h, jdt)
    rng, az, el = razel(rSEZ)
    overpasses = get_overpasses(el, az, rng, dt_ary, rSEZ, jdt=None, rsiteECI=None, rsatECI=None, loc=loc, sat=sat)
    return overpasses


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
        vPEF = np.dot(R, vTEME) + cross(angular_velocity, rPEF)
    else:
        rPEF = np.einsum("ij...,j...->i...", R, rTEME)
        vPEF = (
            np.einsum("ij...,j...->i...", R, vTEME)
            + cross(angular_velocity, rPEF, 0, 0).T
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


if __name__ == "__main__":
    import pickle
    import matplotlib.pyplot as plt
