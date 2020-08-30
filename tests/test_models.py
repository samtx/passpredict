# test_models.py
import pytest

from passpredict import models


def test_find_overpasses_visible_only(init_prediction):
    """
    jd, location, sun, sat are created from fixture
    """
    jd, location, sun, sat = init_prediction
    rho = models.RhoVector(jd, sat, location, sun)
    overpasses = rho.find_overpasses(visible_only=True)
    for overpass in overpasses:
        assert overpass.type == models.PassType.visible
    


# def test_Time_Model():
#     jd0 = julian_date(2020, 7, 1, 12, 34, 56)
#     jd = jd0 + np.linspace(0, 14, 100)
#     t = models.Time(jd=jd)
#     t2 = t[::2]
#     assert_array_equal(t.jd, jd)
#     assert_array_equal(t2.jd, jd[::2])


if __name__ == "__main__":
    pytest.main(['-v', __file__])