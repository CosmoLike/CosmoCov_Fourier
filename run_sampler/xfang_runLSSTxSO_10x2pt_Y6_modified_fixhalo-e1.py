#!/usr/bin/env python

## Photo-z and shear calibration marginalization already absorbed in Covariance

import sys, os
sys.path.append(os.path.dirname(sys.path[0]))

from cosmolike_libs_LSSTxSO_10x2pt_e1 import * 
from schwimmbad import MPIPool

inv=['invcov_Y1_10x2pt','invcov_lsstxso_y6_modified_10x2pt']

data=['10x2pt_LSSTxSO_Y1','10x2pt_fid-e1-0.1']

mask=['...','10x2pt_LSSTxSO_Y6_mask.txt']
# bary=['LPC_6x2pt_LSSTxSO_Y1','LPC_6x2pt_LSSTxSO_Y6']

source_z=['src_LSSTY1','src_LSSTY6'] 

lens_z=['lens_LSSTY1','lens_LSSTY6']

shear_prior=[0.013,0.003] 
delta_z_prior_shear=[0.002,0.001]
delta_z_prior_clustering=[0.002,0.001] ## Not SRD, assumed to be no worse than shear
sigma_z_shear=[0.05,0.05]
sigma_z_clustering=[0.03,0.03]
sigma_z_prior_shear=[0.006,0.003]
sigma_z_prior_clustering=[0.006,0.003] ## Not SRD, assumed to be no worse than shear

nsource_table=[10.7,22.5]  
nlens_table=[13.1,26.8]
area_table=[12300.0,16500.0]

survey_designation=["LSSTxSO_Y1","LSSTxSO_Y6"]
tomo_binning_source=["source_std","source_std"]
tomo_binning_lens=["LSST_gold","LSST_gold"]

model=1
file_source_z = os.path.join(dirname, "zdistris/",source_z[model])
file_lens_z = os.path.join(dirname, "zdistris/",lens_z[model])
data_file = os.path.join(dirname, "datav/",data[model])
cov_file = os.path.join(dirname, "cov/",inv[model])
mask_file = os.path.join(dirname, "datav/",mask[model])
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv[model])
chain_file = os.path.join(dirname, "chains/LSSTxSO_10x2pt_model_%d_modified_fixhalo-e1-0.1" %model)
# bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halomodel".encode('utf-8'))
# initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initbins(25,20.0,7979.0,7979.0,21.0,10,10)

initpriors(shear_prior[model],sigma_z_shear[model],delta_z_prior_shear[model],sigma_z_prior_shear[model],sigma_z_clustering[model],delta_z_prior_clustering[model],sigma_z_prior_clustering[model])
initsurvey(survey_designation[model].encode('utf-8'),nsource_table[model],nlens_table[model],area_table[model])
initgalaxies(file_source_z.encode('utf-8'),file_lens_z.encode('utf-8'),"gaussian".encode('utf-8'),"gaussian".encode('utf-8'),tomo_binning_source[model].encode('utf-8'),tomo_binning_lens[model].encode('utf-8'))
initia("NLA_z".encode('utf-8'),"none".encode('utf-8'))

initfb(1)
# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("10x2pt".encode('utf-8'))
initdatainv(cov_file.encode('utf-8'),data_file.encode('utf-8'),mask_file.encode('utf-8'))
initcmb("so_Y5".encode('utf-8'))

# use modified covariance and skip shearcalib and photo-z sampling
skip_shearcalib_phz_sampling()

#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_10x2_fixhalo(get_N_tomo_shear(),get_N_tomo_clustering(), MG=False, w0wa=False, cov_modified=True)

Nwalker = int(sys.argv[1])
sample_main(sample_params,sigma_z_shear[model],sigma_z_clustering[model],1330,Nwalker,1,chain_file, blind=False, pool=MPIPool())

