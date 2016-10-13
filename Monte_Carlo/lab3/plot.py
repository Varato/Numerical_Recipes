import numpy as np
import pprint
from matplotlib import pyplot as plt
from scipy import interpolate

def analyse_file(file_name, color,label_str):
	file = open(file_name, "r")
	T = []
	C_ave = []
	error = []
	lines = file.readlines()
	for line in lines:
		data = line.strip().split(",")
		T.append(eval(data[0]))
		C_ave.append(1*eval(data[1]))
		error.append(eval(data[2]))
	file.close()
	tck = interpolate.splrep(T, C_ave)
	Tnew = np.arange(T[0],T[-1],0.0001)
	Cnew = interpolate.splev(Tnew, tck)
	C_max = max(Cnew)
	index = list(Cnew).index(C_max)
	Tc = Tnew[index]
	plt.plot(Tnew, Cnew, color+"-")
	plt.errorbar(T, C_ave, yerr=error, fmt=color+"o",label=label_str)
	plt.plot([Tc, Tc], [0, C_max], "k-.")
	return Tc


plt.figure(figsize=[12,9])
plt.title("Honeycomb Ising Model")
plt.axis([0,3.5,0,2000])	
Tc10=analyse_file("result_ave1024.txt","c", "N=1000")
Tc9=analyse_file("result_ave900.txt","r", "N=900")
Tc4=analyse_file("result_ave400.txt", "b", "N=400")
Tc1=analyse_file("result_ave100.txt","m","N=100")
print Tc1, Tc4, Tc9, Tc10

plt.legend(loc="upper right")
plt.xlabel("$k_BT/J$")
plt.ylabel("$C/k_B$")
plt.savefig("result.png")
plt.show()