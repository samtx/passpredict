from passpredictor.timefn_ext import jday2datetime as jday2datetime_cython
from passpredictor.timefn import jday2datetime
import cProfile
import numpy as np

jday = 2450383.09722222  # 1996-10-26 14:20:0
jday_ary = np.array([jday])

# Profile python function
cProfile.run('jday2datetime(jday)')
