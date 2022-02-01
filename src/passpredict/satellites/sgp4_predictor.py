from __future__ import annotations
import datetime


from .base import SatellitePredictorBase, LLH


class SGP4Predictor(SatellitePredictorBase):
    """
    Predictor for satellite overpasses. Uses sgp4 for propagation
    """
    pass