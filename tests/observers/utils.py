import datetime

import pytest


def assert_datetime_approx(dt1: datetime.datetime, dt2: datetime.datetime, delta_seconds: float):
    """
    Compare two python datetimes. Assert difference is <= delta_seconds
    """
    diff = (dt1 - dt2).total_seconds()
    assert diff == pytest.approx(0.0, abs=delta_seconds)