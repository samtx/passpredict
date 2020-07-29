
from collections import namedtuple


EOP = namedtuple('EOP', ['xp', 'yp', 'dUTC1'], defaults=[0., 0., 0.])


def eop(jdt):
    """Return earth orientation parameters for julian date
    
    Returns:
        EOP tuple

    """
    return EOP()