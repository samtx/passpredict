import math

import numpy as np

from .constants import DEG2RAD


def rot1(a):
    """Compute Euler angle rotation matrix, first angle

    Params:
        a : float (n)

    Returns:
        float (3, 3, n)

    References:
        Vallado, Eq. 3-15
        skyfield.functions.rot_x
    """
    c = np.cos(a)
    s = np.sin(a)
    zero = a * 0.0
    one = zero + 1.0
    return np.array(((one, zero, zero), (zero, c, -s), (zero, s, c)))


def rot2(a):
    """Compute Euler angle rotation matrix, second angle
    
    Params:
        a : float (n)

    Returns:
        float (3, 3, n)    
    
    References:
        Vallado, Eq. 3-15
        skyfield.functions.rot_y
    """
    c = np.cos(a)
    s = np.sin(a)
    zero = a * 0.0
    one = zero + 1.0
    return np.array(((c, zero, s), (zero, one, zero), (-s, zero, c)))


def rot3(a):
    """Compute Euler angle rotation matrix, third angle

    Params:
        a : float (n)

    Returns:
        float (3, 3, n)
        
    References:
        Vallado, Eq. 3-15
        skyfield.functions.rot_z
    """
    c = np.cos(a)
    s = np.sin(a)
    zero = a * 0.0
    one = a + 1.0
    return np.array(((c, -s, zero), (s, c, zero), (zero, zero, one)))


def mxm(M1, M2):
    """Matrix times matrix: multiply two NxN matrices.
    
    Reference:
        skyfield.functions
    """
    return np.einsum('ij...,jk...->ik...', M1, M2)


def mxmxm(M1, M2, M3):
    """Matrix times matrix times matrix: multiply 3 NxN matrices together.
    
    Reference:
        skyfield.functions
    """
    return np.einsum('ij...,jk...,kl...->il...', M1, M2, M3)


def mxv(mtx, vec):
    """
    Multiply a matrix by (m) vectors

    Args:
        mtx: float (3, 3, n)
        vec: float (3, n)

    Returns:
        float (3, m)

    Reference:
        skyfield/functions.py, line 23
    """
    return np.einsum('ij...,j...->i...', mtx, vec)


def ecef2sez(r: np.ndarray, phi: float, lmda: float) -> np.ndarray:
    """
    Rotate r vector from ECEF frame to SEZ frame
    Example uses USAF academy as the station
       lmda = -104.0 deg longitude
       phi  =   39.0 deg latitude
    """
    phi_rad = phi * DEG2RAD
    lmda_rad = lmda * DEG2RAD
    ang1 = (90 - phi) * DEG2RAD
    cosang1 = math.cos(ang1)
    sinang1 = math.sin(ang1)
    cosang2 = math.cos(lmda_rad)
    sinang2 = math.sin(lmda_rad)
    M = np.array(
        (
            (cosang1 * cosang2, cosang1 * sinang2, -sinang1),
            (-sinang2, cosang2, 0),
            (sinang1 * cosang2, sinang1 * sinang2, cosang1)
        ),
        dtype=np.float64
    )
    return mxv(M, r)