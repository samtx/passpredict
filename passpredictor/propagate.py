
import numpy as np
from passpredictor.sgp4 import cpropagation

cython_installed = False
try:
    from sgp4.cpropagation import sgp4
except ImportError:
    from sgp4.propagation import sgp4

