import numpy as np
import sys

def gauss(M):
	n = M.shape[0]
	return np.zeros(n)

arq = open(sys.argv[1])
lines = arq.readlines()

A = np.array([[float(x) for x in line.split(' ')] for line in lines[1:]])
print(A)

arq.close()
