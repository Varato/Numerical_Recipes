import numpy as np
import matplotlib.pyplot as plt

x0=-5.
p0=0.5
dt = 0.1
sigma = 1
a = 5*np.sqrt(2)

x_step=0.0001
x=np.arange(-20,20,x_step)
N=len(x)
print N
exit()
global wave_packet
V_exp_factor = []
T_exp_factor = []
wave_packet = []

def Gaussian_wp(x,x0,p0):
	return np.exp(-(x-x0)**2/(2*sigma**2))*np.exp(1j*p0*(x-x0))

def V(x):
	if 0<=x<=a:
		return 1
	else:
		return 0
def normalize():
	global wave_packet
	ww=[]
	for i in range(N):
		ww.append(wave_packet[i].conj()*wave_packet[i])
	prod = np.trapz(ww, x)
	wave_packet = wave_packet/prod


for xx in x:
	wave_packet.append(Gaussian_wp(xx, x0, p0))
normalize()

for xx in x:
	V_exp_factor.append(np.exp(-1j*V(xx)*dt))

for i in range(N):
	if i<N/2:
		p=i*2*np.pi/x_step/N
	else:
		p=(i-N)*2*np.pi/x_step/N
	T_exp_factor.append(np.exp(-1j*p**2*dt/2))


def evolve(n):
	global wave_packet
	for t in range(n):
		for i in range(N):
			wave_packet[i] *= V_exp_factor[i]
		p_wave_packet = np.fft.fft(wave_packet)/np.sqrt(2*np.pi)
		for i in range(N):
			p_wave_packet[i] *= T_exp_factor[i]
		wave_packet = np.fft.ifft(p_wave_packet)/np.sqrt(2*np.pi)
		normalize()


# plt.figure(figsize=[19.24,10.80])
# plt.subplot(421)
# # plt.plot([0,0,a,a,0],[-1,1,1,-1], 'k-')
# plt.plot(x, wave_packet)

# evolve(15)
# plt.subplot(422)
# plt.plot(x, wave_packet)

# evolve(20)
# plt.subplot(423)
# # plt.plot([0,0,a,a],[-1,1,1,-1], 'k-')
# plt.plot(x, wave_packet)

# evolve(20)
# plt.subplot(424)
# # plt.plot([0,0,a,a],[-1,1,1,-1], 'k-')
# plt.plot(x, wave_packet)

# evolve(20)
# plt.subplot(425)
# # plt.plot([0,0,a,a],[-1,1,1,-1], 'k-')
# plt.plot(x, wave_packet)


evolve(20)
plt.subplot(426)
# plt.plot([0,0,a,a],[-1,1,1,-1], 'k-')
plt.plot(x, wave_packet)

plt.show()
