import numpy as np

# Constants
R_EARTH  = np.float64(6378.137) # [km] Mean equatorial radius
R2_EARTH = np.longdouble(40680631.59076899)  # [km2] mean equatorial radius squared
e_EARTH  = np.longdouble(0.081816221456)  # Earth eccentricity
e2_EARTH = np.longdouble(0.006694385000)  # Earth eccentricity squared
MU = np.longdouble(398600.4418)  # [km3/(solar s)2] gravitational parameter
J2 = np.longdouble(0.0010826267)


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

def nu2anomaly(nu, e):
    """
    Vallado, p.77
    """
    nu_rad = np.radians(nu)
    if e < 1.0:
        # sinE = np.sin(nu_rad)*np.sqrt(1-e**2)/(1+e*np.cos(nu_rad))
        # E = np.arcsin(sinE)
        cosE = (e + np.cos(nu_rad))/(1 + e*np.cos(nu_rad))
        E = np.arccos(cosE)
    return np.degrees(E)

def anomaly2nu(e, E):
    """
    Vallado, alg 6, p.77
    """
    E_rad = np.radians(E)
    if e < 1.0:
        cosNu = (np.cos(E_rad) - e)/(1 - e*np.cos(E_rad))
        nu = np.arccos(cosNu)
    return np.degrees(nu)

def coe2rv():
    pass

def kepEqtnE(M, e):
    """
    Vallado, Alg 2, p.65
    """
    M_rad = np.radians(M)
    if (-np.pi < M_rad < 0) or (M_rad > np.pi):
        E0 = M_rad - e
    else:
        E0 = M_rad + e
    atol = 1e-8
    E1 = 1.
    k = 0
    while k < 1000:
        E1 = E0 + (M_rad - E0 + e*np.sin(E0))/(1 - e*np.cos(E0)) 
        if np.abs(E1 - E0) < atol:
            break
        E0 = E1
        k += 1
    return np.degrees(E1)

def pkepler_coe(a0, e0, i0, Omega0, w0, nu0, u, lmda_true, what_true, dt, ndt=0., nddt=0.):
    """
    Vallado, Alg 65, p.691

    Propagate object from orbital elements
    """

    if e != 0.0:
        E0 = nu2anomaly(nu0, e0)

    M0 = E0 - e0*np.sin(np.radians(E0))
    p0 = a0*(1 - e0**2)
    n0 = np.sqrt(MU/a0**3)

    # Update for permutations
    a = a0 - 2*a0/(3*n0)*ndt*dt
    e = e0 - 2*(1 - e0)/(3*n0)*ndt*dt
    Omega = Omega0 - 3*n0*R2_EARTH*J2/(2*p0**2)*np.cos(np.radians(i0))*dt
    w = w0 + 3*n0*R2_EARTH*J2/(4*p0**2)*(4 - 5*np.sin(np.radians(i0))**2)*dt
    M = M0 + n0*dt + ndt/2.*dt**2 + nddt/6.*dt**3
    p = a*(1-e**2)

    E = kepEqtnE(M, e)
    if e != 0.0:
        nu = anomaly2nu(E, e)
    else:
        u, lmda_true = E, E


def sgp4():
    pass