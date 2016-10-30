import numpy as np
N=8
h=np.pi/4
xx=np.zeros(N)
x=[]
for i in range(N):
	xx[i]=i*h	
	x.append(np.exp(1j*xx[i]))
print x
y = np.fft.fft(x).round()
print y
print np.fft.ifft(y).round()
