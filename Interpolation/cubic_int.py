import matplotlib.pyplot as plt
import numpy as np

x = [-1, 1, 2, 4]
y = [1.25, 2, 3, 0]

plt.plot(x, y, 'ro', markersize=5)
plt.axis([-2,5,-1,4])

c = [[0.089, 0.268, 0.286, 1.357], [-0.625, 2.411, -1.857, 2.071], [0.223, -2.679, 8.321, -4.714]]

def S(x, i):
	cc = c[i-1]
	s = cc[0]*x**3 + cc[1]*x**2 + cc[2]*x + cc[3]
	return s
x1 = np.linspace(-1, 1, 50)
x2 = np.linspace(1, 2, 50)
x3 = np.linspace(2, 4, 50)
plt.plot(x1, S(x1, 1), "b-")
plt.plot(x2, S(x2, 2), "b-")
plt.plot(x3, S(x3, 3), "b-")
plt.show()
