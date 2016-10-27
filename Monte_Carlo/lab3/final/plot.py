import numpy as np
import pprint
from matplotlib import pyplot as plt
from scipy import interpolate

def analyse_file(file_name, color,label_str):
	N = eval(label_str.split("=")[1])
	file = open(file_name, "r")
	T = []
	C_ave = []
	error = []
	lines = file.readlines()
	for line in lines:
		data = line.strip().split(",")
		T.append(eval(data[0]))
		C_ave.append(1.*eval(data[1])/N)
	#	error.append(eval(data[2]))
	file.close()
	tck = interpolate.splrep(T, C_ave)
	Tnew = np.arange(T[0],T[-1],0.0001)
	Cnew = interpolate.splev(Tnew, tck)
	C_max = max(Cnew)
	index = list(Cnew).index(C_max)
	Tc = Tnew[index]
	plt.plot(Tnew, Cnew, color+"-")
	plt.plot(T, C_ave, color+"o", label=label_str)
	#plt.errorbar(T, C_ave, yerr=error, fmt=color+"o",label=label_str)
	plt.plot([Tc, Tc], [0, C_max], "k-.")
	return Tc

def energy_plot(file_name, color,label_str):
	file = open(file_name, "r")
	T = []
	E_ave = []
	error = []
	lines = file.readlines()
	for line in lines:
		data = line.strip().split(",")
		T.append(eval(data[0]))
		E_ave.append(1*eval(data[3]))
		error.append(eval(data[4]))
	file.close()
	plt.plot(T, E_ave, color+"-.")
	plt.errorbar(T, E_ave, yerr=error, fmt=color+"o",label=label_str)

plt.figure(1,figsize=[15,10])
# plt.subplot(121)
# Tc10=analyse_file("result1024.txt","c", "N=1024")
Tc9=analyse_file("result900.txt","r", "N=900")
Tc4=analyse_file("result400.txt", "b", "N=400")
Tc1=analyse_file("result100.txt","m","N=100")
print Tc1, Tc4, Tc9#, Tc10
plt.legend(loc="upper right")
plt.xlabel("$k_BT/J$")
plt.ylabel("$c/k_B$")
plt.title("Specific Capacity Plot")
plt.savefig("result1.png")

plt.figure(2,figsize=[15,10])
# energy_plot("result1024.txt","c", "N=1024")
energy_plot("result900.txt","r", "N=900")
energy_plot("result400.txt","b", "N=400")
energy_plot("result100.txt","m", "N=100")
plt.legend(loc="lower right")
plt.xlabel("$k_BT/J$")
plt.ylabel("$<E>/J$")
plt.title("Energy Plot")
plt.savefig("result2.png")
plt.show()
