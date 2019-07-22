import numpy as np

# Constants
R_EARTH  = np.float64(6378.137) # [km] Mean equatorial radius
e_EARTH  = np.float64(0.081816221456)  # Earth eccentricity
e2_EARTH = np.float64(0.006694385000)  # Earth eccentricity squared


def coe2rv():
    pass

def site_declination_and_K(phi_gd, h_ellp):
    """
    Vallado, Eq.3-7 
    Get declination and K vectors for site on Earth
    Note: currently only precise to 0.1 km
    """
    # Eq.(3-7)
    phi_gd_rad = np.radians(np.float64(phi_gd))  # convert [deg] to [rad]
    h_ellp_km = np.float64(h_ellp)*(1/1000)      # convert [m] to [km]
    C = R_EARTH / np.sqrt(1 - e2_EARTH*np.sin(phi_gd_rad)**2)
    S = C*(1-e2_EARTH)
    r_delta = (C + h_ellp_km)*np.cos(phi_gd_rad)
    r_K = (S + h_ellp_km)*np.sin(phi_gd_rad)
    return (r_delta, r_K)


def site_ECEF(phi_gd, lmda, h_ellp):
    """
    Algorithm 51 from Vallado, p.430

    """

    r_delta, r_K = site_declination_and_K(phi_gd, h_ellp)

    lmda_rad = np.radians(lmda)

    r_site_ECEF = np.array([
        r_delta*np.cos(lmda_rad),
        r_delta*np.sin(lmda_rad),
        r_K
    ])

    return r_site_ECEF



def fk5():
    pass

def pkepler():
    pass

def sgp4():
    pass