import numpy as np
import matplotlib.pyplot as plt

wave_packet=np.loadtxt("result")
x=np.arange(-20,20,0.001)
for i in range(21):
	plt.figure(i)
	plt.plot(x, wave_packet[i])
	plt.savefig("p%d"%i)
