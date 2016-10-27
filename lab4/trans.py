import numpy as np
import matplotlib.pyplot as plt

# x=np.linspace(-5,5,50)
# plt.plot(x, np.sinh(x))
# plt.show()
V0 = 1.
a = np.sqrt(2)*5.

def T(E):
	k = np.sqrt(2*E)
	kapa = np.sqrt(2*np.abs(V0-E)) 
	up = 4*k**2*kapa**2
	if E<=1.:
		down = 4*k**2*kapa**2*np.cosh(kapa*a)**2 + (kapa**2 - k**2)**2*np.sinh(kapa*a)**2
	else:
		down = 4*k**2*kapa**2*np.cos(kapa*a)**2 + (kapa**2 + k**2)**2*np.sin(kapa*a)**2
	return up/down

E=np.linspace(0,5,10000)
TT=[]
for ee in E:
	TT.append(T(ee))
plt.figure(figsize=[14,10])
plt.plot(E,TT,linewidth=2)
plt.plot([2,2],[0,T(2.)], 'k-.')
plt.title("Transmission Probability")
plt.xlabel("$E/V_0$")
plt.ylabel("$T$")
plt.savefig("transprob.png")
plt.show()
