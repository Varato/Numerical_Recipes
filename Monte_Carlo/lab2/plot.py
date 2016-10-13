import numpy
from matplotlib import pyplot as plt
import pprint

file_name = "result.txt"
file = open(file_name, "r")
rAB = []
E_ave = []
Error = []
a0 = 0.529#(17721067) Angstrom
e_d_a0 = 27.2 #eV
for i in range(30):
	line = file.readline()
	data = line.strip().split(",")
	rAB.append(a0*eval(data[0]))
	E_ave.append(e_d_a0*eval(data[1]))
	Error.append(e_d_a0*eval(data[2]))
file.close()
plt.figure(figsize=[15,10])
plt.errorbar(rAB, E_ave, yerr=Error)
plt.plot(rAB, E_ave, "r-")
plt.xlabel("$r_{AB} / \AA$")
plt.ylabel("$E / eV$")
plt.savefig("result.png")
plt.show()
