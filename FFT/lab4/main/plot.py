import numpy as np
import matplotlib.pyplot as plt

wave_packet=np.loadtxt("result")
x=np.arange(-20,20,0.001)
l=wave_packet.shape[0]
for i in range(l):
	plt.figure()
	plt.plot(x, wave_packet[i])
	plt.savefig("p%d"%i)
