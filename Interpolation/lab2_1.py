from ctypes import *
import os
import numpy as np
import matplotlib.pyplot as plt

mod=cdll.LoadLibrary(os.getcwd()+'/polint_mod.so')
mod.interpolate.restype=c_float
xa = [-1., 1., 2., 4.]
ya = [1.25, 2., 3., 0]
x = np.linspace(-1,4,80)
y = []

for xx in x:
	y.append(mod.interpolate(c_float(xx)))
plt.plot(x,y, color="cornflowerblue",label="interpolation")
plt.plot(xa,ya,'ro', markersize=8, label="data points")
plt.legend(loc="upper, right")
plt.axis([-1.5,4.5, -0.5,3.5])
plt.grid(True)
plt.show()


