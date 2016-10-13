#/Users/xinchen/anaconda/bin/python
import numpy as np
from numpy import linalg as la
import pprint

A = np.array([[6,2],[2,2]])

U,S,V = la.svd(A)

a=np.array([[4,4],[9,9]])
print np.sqrt(a)

