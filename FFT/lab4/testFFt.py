import numpy as np
import matplotlib.pyplot as plt

N=8
xx=[]
kk=[]
tmp=0
tmp1=0
for i in range(N):
	xx.append(tmp)
	kk.append(tmp1)
	tmp += np.pi/4
	tmp1 += 2*np.pi*4/np.pi/N
x = np.sin(xx)
print sum(x)
k = np.fft.fft(x)
plt.subplot(121)
plt.plot(kk,k)
plt.subplot(122)
plt.plot(xx,np.fft.ifft(k))
plt.show()