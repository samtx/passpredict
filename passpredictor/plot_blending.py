# test.py

import numpy as np
import matplotlib.pyplot as plt
from passpredictor import parabolic_blending, compute_alphas

C = lambda x: -4 * x ** 3 + 3 * x ** 2 + 25 * x + 6

x = np.linspace(-5, 5)
y = C(x)
t0 = -5
tf = 5
t = np.linspace(t0, tf, 4)
yt = C(t)

trise, tset = parabolic_blending(t, yt)
n, m = len(trise), len(tset)
troots = np.empty(n + m)
troots[0:n] = trise[:]
troots[n : n + m] = tset[:]
ytroots = C(troots)
plt.plot(x, y, ":r", label="C(x)")
plt.scatter(t, yt, facecolors=None, label="C(t)")
plt.scatter(troots, ytroots, facecolors=None, label="C(troots)")
plt.legend()
plt.grid()
plt.show()
