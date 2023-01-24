import numpy as np
import sys
import os

if __name__ == '__main__':
	filename = sys.argv[1]
	maskname = os.path.splitext(filename)[0] + '_mask.txt'

	dv = np.loadtxt(filename)[:,1]
	mask = np.ones(dv.size, dtype=int)
	mask[dv==0.0] = 0

	np.savetxt(maskname, np.c_[np.arange(dv.size), mask], fmt="%d %d")
	print("mask file %s saved!"%(maskname))
