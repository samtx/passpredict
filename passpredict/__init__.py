try:
    from .predict import Location, Point
except ImportError:
    raise ImportError('Must compile with Cython first')

from .tle import OMM