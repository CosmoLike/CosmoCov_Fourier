import numpy as np 

infile = 'cmb_lmax3000.txt'

outfile = 'cmb_lmax3000_new.txt'

data = np.loadtxt(infile)
ell = data[:,0]
noise = data[:,1]

noise *= ell*(ell+1.)/4.

np.savetxt(outfile, np.c_[ell, noise], fmt='%.15e %.15e')