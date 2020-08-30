# test_models.py
import pytest
import numpy as np

from passpredict import models


def test_find_overpasses_visible_only(init_find_overpasses):
    """
    jd, location, sun, sat are created from fixture
    """
    jd, location, sun, sat = init_find_overpasses
    rho = models.RhoVector(jd, sat, location, sun)
    overpasses = rho.find_overpasses(visible_only=True)
    for overpass in overpasses:
        assert overpass.type == models.PassType.visible
    


def test_SatPredictModel_slice():
    n = 40
    sat = models.SatPredictData(
        id=25544,
        rECEF=np.linspace(0, n*3, n*3).reshape((3,n)),
        illuminated=np.ones(n, dtype=bool)
    )
    slc = slice(10)
    sat2 = sat[slc]
    assert sat2.id == sat.id
    assert np.array_equal(sat2.rECEF, sat.rECEF[:,slc])
    assert np.array_equal(sat2.illuminated, sat.illuminated[slc])


def test_SunPredictModel_slice():
    n = 40
    sun = models.SunPredictData(
        rECEF=np.linspace(0, n*3, n*3).reshape((3,n)),
    )
    slc = slice(10)
    sun2 = sun[slc]
    assert np.array_equal(sun2.rECEF, sun.rECEF[:,slc])
    

def test_SunPredictModel_slice_array_view():
    n = 40
    sun = models.SunPredictData(
        rECEF=np.linspace(0, n*3, n*3).reshape((3,n)),
    )
    slc = slice(10)
    sun2 = sun[slc]
    assert np.may_share_memory(sun2.rECEF, sun.rECEF)
    assert np.may_share_memory(sun2.rECEF, sun.rECEF[:, slc])


def test_SatPredictModel_slice_array_view():
    n = 40
    sat = models.SatPredictData(
        id=25544,
        rECEF=np.linspace(0, n*3, n*3).reshape((3,n)),
        illuminated=np.ones(n, dtype=bool)
    )
    slc = slice(10)
    sat2 = sat[slc]
    assert id(sat2.id) == id(sat.id)
    assert np.may_share_memory(sat2.rECEF, sat.rECEF[:, slc])
    assert np.may_share_memory(sat2.rECEF, sat.rECEF)
    assert np.may_share_memory(sat2.illuminated, sat.illuminated[slc])
    assert np.may_share_memory(sat2.illuminated, sat.illuminated)

# def test_Time_Model():
#     jd0 = julian_date(2020, 7, 1, 12, 34, 56)
#     jd = jd0 + np.linspace(0, 14, 100)
#     t = models.Time(jd=jd)
#     t2 = t[::2]
#     assert_array_equal(t.jd, jd)
#     assert_array_equal(t2.jd, jd[::2])


if __name__ == "__main__":
    pytest.main(['-v', __file__])