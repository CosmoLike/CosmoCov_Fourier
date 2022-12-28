#!/usr/bin/env python

## Photo-z and shear calibration marginalization already absorbed in Covariance

import sys, os
sys.path.append(os.path.dirname(sys.path[0]))

from cosmolike_libs_RomanxSO_10x2pt import * 
from schwimmbad import MPIPool

survey_designation="RomanxSO"
probes = ['10x2pt', '3x2pt', '6x2pt']
i = int(sys.argv[1])
probe = probes[i]

inv=f'invcov_romanxso_modified_{probe}'

data=f'{probe}_{survey_designation}'

mask=f'{probe}_{survey_designation}_mask.txt'

source_z='zdistri_WFIRST_LSST_lensing_fine_bin_norm'

lens_z='zdistri_WFIRST_LSST_clustering_fine_bin_norm'

sigma_z_shear=0.01
sigma_z_clustering=0.01

lmax_shear = 4000.

# No need:
shear_prior=0.003
delta_z_prior_shear=0.001
delta_z_prior_clustering=0.001 ## Not SRD, assumed to be no worse than shear
sigma_z_prior_shear=0.003
sigma_z_prior_clustering=0.003 ## Not SRD, assumed to be no worse than shear
nsource_table=66.0  
nlens_table=51.0
area_table=2000.0
####

tomo_binning_source="source_std"
tomo_binning_lens="WF_SN10"

file_source_z = os.path.join(dirname, "zdistris/",source_z)
file_lens_z = os.path.join(dirname, "zdistris/",lens_z)
data_file = os.path.join(dirname, "datav/",data)
cov_file = os.path.join(dirname, "cov/",inv)
mask_file = os.path.join(dirname, "datav/",mask)
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv)
chain_file = os.path.join(dirname, f"chains/{survey_designation}_{probe}_modified")
# bary_file=os.path.join(dirname, "baryons/",bary)

initcosmo("halomodel".encode('utf-8'))
# initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initbins(25,20.0,7979.0,lmax_shear,21.0,10,10)

initpriors(shear_prior,sigma_z_shear,delta_z_prior_shear,sigma_z_prior_shear,sigma_z_clustering,delta_z_prior_clustering,sigma_z_prior_clustering)
initsurvey(survey_designation.encode('utf-8'),nsource_table,nlens_table,area_table)
initgalaxies(file_source_z.encode('utf-8'),file_lens_z.encode('utf-8'),"gaussian".encode('utf-8'),"gaussian".encode('utf-8'),tomo_binning_source.encode('utf-8'),tomo_binning_lens.encode('utf-8'))
initia("NLA_z".encode('utf-8'),"none".encode('utf-8'))

initfb(1)
# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes(probe.encode('utf-8'))
initdatainv(cov_file.encode('utf-8'),data_file.encode('utf-8'),mask_file.encode('utf-8'))
initcmb("so_Y5".encode('utf-8'))

# use modified covariance and skip shearcalib and photo-z sampling
skip_shearcalib_phz_sampling()

#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_10x2_allsys(get_N_tomo_shear(),get_N_tomo_clustering(), MG=False, w0wa=False, cov_modified=True)

Nwalker = int(sys.argv[2])
sample_main(sample_params,sigma_z_shear,sigma_z_clustering,8000,Nwalker,1,chain_file, blind=False, pool=MPIPool())
