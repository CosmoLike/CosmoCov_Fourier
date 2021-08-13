import numpy as np 
import pandas as pd 

filenames = ["Apr17_mv_nlkk_deproj0_BASELINE_fsky_16000_iterOn.csv", "Apr17_mv_nlkk_deproj0_GOAL_fsky_16000_iterOn.csv"]
savenames = ["so_baseline_nlkk_lmax3000.txt", "so_gold_nlkk_lmax3000.txt"]

for i in range(2):
	d = pd.read_csv(filenames[i], delim_whitespace=True)

	ell = np.array(d.ix[:,0]).astype(int)
	nlkk = np.array(d.ix[:,1])
	print(ell.size)
	print(nlkk.size)

	np.savetxt(savenames[i], np.c_[ell, nlkk])