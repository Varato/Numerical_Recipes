import numpy as np
import matplotlib.pyplot as plt

def wave_packet(x,x0,p0):
	return np.exp(-(x-x0)**2/2)*np.exp(1j*p0*(x-x0))

x0=0.
p0=10.
x=np.linspace(-10,10,1000)
f=[]
for xx in x:
	f.append(wave_packet(xx, x0, p0).real)

plt.plot(x,f)
plt.show()