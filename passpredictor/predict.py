import numpy as np
from numpy import dot, cross
from numpy.linalg import norm
import math
from math import sqrt, sin, cos, cosh, acosh, tan, atan, acos, radians, degrees, pi
import datetime
from passpredictor.rotations import rot1, rot2, rot3, theta_GMST1982
from passpredictor.solar import sun_pos, is_sat_illuminated
from passpredictor.constants import (
    R_EARTH, R2_EARTH, e_EARTH, e2_EARTH, MU, J2, J2000, AU_M, AU_KM, ASEC360,
    DAY_S, ASEC2RAD, DEG2RAD, RAD2DEG, tau
)


def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    phi_gd_rad = radians(np.float64(phi_gd))  # convert [deg] to [rad]
    h_ellp_km = np.float64(h_ellp) * (1 / 1000)  # convert [m] to [km]
    C = R_EARTH / np.sqrt(1 - e2_EARTH * sin(phi_gd_rad) ** 2)
    S = C * (1 - e2_EARTH)
    r_delta = (C + h_ellp_km) * cos(phi_gd_rad)
    r_K = (S + h_ellp_km) * sin(phi_gd_rad)
    return (r_delta, r_K)


def site_ECEF(phi_gd, lmda, h_ellp):
    """Compute ECEF coordinates for tracking site on Earth

    References:
        Vallado, Algorithm 51, p.430
    """
    r_delta, r_K = site_declination_and_K(phi_gd, h_ellp)
    lmda_rad = radians(lmda)
    r_site_ECEF = np.array([r_delta * cos(lmda_rad), r_delta * sin(lmda_rad), r_K])
    return r_site_ECEF


def site_ECEF2(phi_gd, lmda, h_ellp):
    """Compute ECEF coordinates for tracking site on Earth

    Args:
        phi_gd: (float) geodetic latitutde of site in degrees
        lmda: (float) east longitude of site in degrees
        h_ellp: (float) height above the reference ellipse in meters

    References:
        Vallado, p. 428, Eq. 7-1
    """
    phi_gd_rad = phi_gd * DEG2RAD
    lmda_rad = lmda * DEG2RAD
    cosphi = math.cos(phi_gd_rad)
    sinphi = math.sin(phi_gd_rad)
    C = R_EARTH / math.sqrt(1 - e2_EARTH * (sinphi ** 2))
    S = C * (1 - e2_EARTH)
    h_ellp *= 0.001  # convert to km
    tmp = (C + h_ellp) * cosphi
    r_site_ECEF = np.array(
        [tmp * math.cos(lmda_rad), tmp * math.sin(lmda_rad), (S + h_ellp) * sinphi]
    )
    return r_site_ECEF


def rng_el(r):
    """Get range and elevation from SEZ vector"""
    rng = np.linalg.norm(r, axis=1)
    el = np.arcsin(r[2] / rng)
    el *= RAD2DEG
    return el, rng


def razel(r):
    """Get range, azimuth, and elevation from SEZ vector"""
    rng = np.linalg.norm(r, axis=0)
    el = np.arcsin(r[2] / rng) * RAD2DEG
    az = (np.arctan2(r[0], r[1]) + pi * 0.5) * RAD2DEG
    idx = np.all([r[0] < 0,r[1] < 0], axis=0)
    az[idx] %= 360
    return rng, az, el


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


def azm(s, e):
    """Compute azimuth from topocentric horizon coordinates SEZ
    Args:
        s : south vector from SEZ coordinate
        e : east vector from SEZ coordinate
    Output:
        azimuth angle in radians with 0 at north.
    """
    out = np.arctan2(s, e) + pi * 0.5
    if s < 0 and e < 0:
        out = out % (2 * pi)
    return out


def elev(z, rhomag):
    """Compute elevation angle from topocentric horizon coordinates SEZ
    Args:
        z : Z vector from SEZ coordinate
        rhomag : magnitude of SEZ coordinate vector
    Output:
        elevation angle in radians with 0 at horizon, pi/2 straight up
    """
    return np.arcsin(z / rhomag)


def get_overpass_idx(el):
    """Compute overpasses based on elevation angle and return indecies
    Args:
        el : float (n), elevation angle in degrees
    Returns:
        overpasses : list, list of indecies of overpasses
    """
    el0 = el[:-1]
    el1 = el[1:]
    el_change_sign = (el0*el1 < 0)
    # Find the start of an overpass
    start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]
    # Find the end of an overpass
    end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]

    # Iterate over start/end indecies and gather inbetween indecies
    overpasses = np.empty(start_idx.size, dtype=object)
    for j, idx in enumerate(start_idx):
        # Store indecies of overpasses in a list
        overpasses[j] = np.arange(idx, end_idx[j]+1, dtype=int)
    return overpasses


def get_overpasses(el, azm, rng, dt_ary, rhoSEZ, sat, loc):
    el0 = el[:-1]
    el1 = el[1:]
    el_change_sign = (el0*el1 < 0)
    # Find the start of an overpass
    start_idx = np.nonzero(el_change_sign & (el0 < el1))[0]
    # Find the end of an overpass
    end_idx = np.nonzero(el_change_sign & (el0 > el1))[0]
    # print(f'start shape {start_idx.shape}  end shape = start shape {end_idx.shape}')

    # Iterate over start/end indecies and gather inbetween indecies
    overpasses = np.empty(start_idx.size, dtype=object)
    for j in range(start_idx.size):
        # Store indecies of overpasses in a list
        idx0 = start_idx[j]
        idxf = end_idx[j]
        overpass_idx = np.arange(idx0, idxf+1, dtype=int)
        idxmax = np.argmax(el[overpass_idx])
        start_pt = Point(
            dt_ary[idx0],
            az[idx0],
            el[idx0],
            rng[idx0]
        ),
        max_pt = Point(
            dt_ary[idxmax],
            az[idxmax],
            el[idxmax],
            rng[idxmax]
        )
        end_pt = Point(
            dt_ary[idxf],
            az[idxf],
            el[idxf],
            rng[idxf]
        ),
        overpass = Overpass(
            start_pt,
            max_pt,
            end_pt,
            dt_ary[overpass_idx],
            rhoSEZ[:,overpass_idx]
        )
        overpass.location = loc
        overpasses[j] = overpass

    return overpasses


def site_rotations(lat, lon, h, rECEF):
    rsiteECEF = site_ECEF(lat, lon, h)
    rho = rECEF - np.atleast_2d(rsiteECEF).T
    rSEZ = rotations.ecef2sez(rho, lat, lon)
    rng, az, el = razel(rSEZ)
    return rng, az, el


def realtime_compute(lat, lon, h, rECEF, jdt):
    rsiteECEF = predict.site_ECEF(lat, lon, h)
    rho = rECEF - np.atleast_2d(rsiteECEF).T
    rSEZ = rotations.ecef2sez(rho, lat, lon)
    rng, az, el = predict.razel(rSEZ)
    np.savez('out2.npz', az=az, el=el, rng=rng, rSEZ=rSEZ, jdt=jdt)


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
