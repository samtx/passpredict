
# from passpredict.sgp4io import twoline2rv, Satellite, wgs72, wgs84
import numpy as np

# use_cython = False
# try:
#     from passpredict._sgp4 import sgp4 as sgp4_pyx
#     use_cython = True
# except ImportError:
#     pass

# from passpredict.sgp4 import sgp4


ISS_TLE = (
"1 25544U 98067A   19293.90487327  .00016717  00000-0  10270-3 0  9034",
"2 25544  51.6426  97.8977 0006846 170.6875 189.4404 15.50212100 34757"
)

# def test_sgp4_tle():
#     print(ISS_TLE)
#     tle1, tle2 = ISS_TLE
#     print(tle1, tle2)
#     satrec = twoline2rv(tle1, tle2, wgs84)
#     min_per_day = 1440
#     total_days = 1
#     dt_per_min = 10
#     t = np.linspace(0, min_per_day*total_days, min_per_day*total_days*dt_per_min)
#     r = np.empty((t.size, 3))
#     v = np.empty((t.size, 3))
#     for i in range(t.size):
#         if i % 43200 == 0:
#             print(f'i = {i}')
#         ri, vi = sgp4(satrec, t[i], wgs84)
#         r[i] = ri
#         v[i] = vi
#     # print(satrec)
#     assert True
