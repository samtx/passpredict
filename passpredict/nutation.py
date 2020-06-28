# nutation.py
from dataclasses import dataclass

import numpy as np

from .constants import DEG2RAD, ASEC2RAD


@dataclass
class Nutation1980Params:
    eclip: float = None
    M_moon: float = None
    M_sun: float = None
    u_M_moon: float = None
    D_sun: float = None
    omega_moon: float = None


# Largest IAU 1980 Nutation Coefficients
# Reference: Vallado book, p. 1043, Table D-6
nut80_i = np.array([1, 9, 31, 2, 10, 32, 11, 33, 34, 12, 35, 13, 36, 38, 37], dtype=np.int64)
nut80_A = np.array([-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, -158, 129, 123, 63, 63], dtype=np.float64)
nut80_B = np.array([-174.2,-1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0.0, -0.5, 0.0, 0.1, 0.0, 0.1, 0.0], dtype=np.float64)
nut80_C = np.array([92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, -1, -70, -53, -33, -2], dtype=np.int64)
nut80_D = np.array([8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
nut80_an = np.array(
    [
        [0, 0, 0, 0, 1],
        [0, 0, 2, -2, 2],
        [0, 0, 2, 0, 2],
        [0, 0, 0, 0, 2],
        [0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 2, -2, 2],
        [0, 0, 2, 0, 1],
    ],
    dtype=np.int64
)


def nut80_fundamental_arguments(tt):
    """IAU 1980 fundamental arguments, Delaunay parameters

    Reference: 
        Vallado book, Eq. 3-82, p. 225
        McCarthy, 1992, p.23

    """
    r = 1296000.0  # arcseconds in 360 degrees
    # r = 360
    print(ASEC2RAD)
    ASEC2DEG = 1/3600.0

    
    ERFA_DAS2R = 4.848136811095359935899141e-6  # Arcseconds to radians   
    ERFA_TURNAS = 1296000.0  # Arcseconds in a full circle 
    ERFA_D2PI = 6.283185307179586476925287 # two pi

    # void eraAnpm() normalizes angle to range -pi <= a <= +pi

    # Mean anomaly of the Moon (l)
    M_moon = (134.96340251 + (((-0.00024470*tt + 0.051635)*tt + 31.8792)*tt + 1717915923.2178)*tt*ASEC2DEG) % 360.000000
    
    ttt = tt
    M_moon = (((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) / 3600.0 + 134.96298139) % 360.
    M_sun = (((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333) % 360.
    u_M_moon = (((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028) % 360.
    D_sun = (((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306) % 360.
    omega_moon =  (((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222) % 360.
    # # M_moon = (134.96340251 + ((1.4343e-5*tt + 0.0088553)*tt + 198.8675605 + 1325*r)*tt*ASEC2DEG) % 360.0
    # M_sun = (357.52910918 + ((3.8e-8*tt - 0.0001537)*tt + 359.0502911 + 99*r)*tt*ASEC2DEG) % 360.0
    # u_M_moon = (93.27209062 + ((-2.88e-7*tt - 0.0035420)*tt + 82.0174577 + 1342*r*ASEC2RAD)*tt) % 360.0
    # D_sun = (297.85019547 + ((1.831e-6*tt - 0.0017696)*tt + 307.1114469 + 1236*r*ASEC2RAD)*tt) % 360.0
    # omega_moon = (125.04455501 + ((2.139e-6*tt + 0.0020756)*tt - 134.1361851 - 5*r*ASEC2RAD)*tt) % 360.0
    return Nutation1980Params(
        M_moon=M_moon,
        M_sun=M_sun,
        u_M_moon=u_M_moon,
        D_sun=D_sun,
        omega_moon=omega_moon
    )


# l = ((((0.064) * ttt + 31.310) * ttt + 171791
# 5922.6330) * ttt) / 3600.0 + 134.9629 8139;
# l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333;
# f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028;
# d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306;
# omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222;