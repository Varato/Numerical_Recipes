import numpy
from matplotlib import pyplot as plt
import pprint

file_name = "result2.txt"
file = open(file_name, "r")
T = []
C = []
E = []
for i in range(15):
	line = file.readline()
	data = line.strip().split(",")
	T.append(eval(data[0]))
	C.append(eval(data[1]))
	E.append(eval(data[2]))
file.close()
plt.figure(figsize=[15,10])
# plt.errorbar(rAB, E_ave, yerr=Error)
plt.plot(T, C, "r-.o")
#plt.plot(T, E, "b-.*")
plt.xlabel("$k_BT/J$")
plt.ylabel("$C/k_B$")
plt.savefig("result.png")
plt.show()
