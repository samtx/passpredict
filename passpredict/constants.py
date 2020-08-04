import numpy as np

# Constants
R_EARTH = 6378.137  # [km] Mean equatorial radius
R2_EARTH = 40680631.59076899  # [km2] mean equat. radius sq.
e_EARTH = 0.081816221456  # Earth eccentricity
e2_EARTH = 0.006694385000  # Earth eccentricity squared
MU = 398600.4418  # [km3/(solar s)2] gravitational parameter
J2 = 0.0010826267
J2000 = 2451545.0
DJC = 36525.0  # days per Julian century

# Various constants required by Skyfield
AU_M = 149597870700  # per IAU 2012 Resolution B2
AU_KM = 149597870.700
ASEC360 = 1296000.0
DAY_S = 86400.0

# Angles.
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321
tau = 6.283185307179586476925287  # lower case, for symmetry with math.pi
