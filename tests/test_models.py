# test_models.py
import datetime

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from passpredict import models
from passpredict.timefn import julian_date


# def test_Time_Model():
#     jd0 = julian_date(2020, 7, 1, 12, 34, 56)
#     jd = jd0 + np.linspace(0, 14, 100)
#     t = models.Time(jd=jd)
#     t2 = t[::2]
#     assert_array_equal(t.jd, jd)
#     assert_array_equal(t2.jd, jd[::2])


if __name__ == "__main__":
    pytest.main(['-v', __file__])