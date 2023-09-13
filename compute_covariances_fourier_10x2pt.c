#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

#include <fftw3.h>


#include "../cosmolike_core/cfftlog/cfftlog.h"
#include "../cosmolike_core/cfftlog/utils.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo_fast.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D_extend.c"
// #include "../cosmolike_core/theory/covariances_fourier.c" // deprecated in this application!!
#include "../cosmolike_core/theory/CMBxLSS_fourier.c"
// #include "../cosmolike_core/theory/covariances_CMBxLSS_fourier.c" // deprecated in this application!!
#include "../cosmolike_core/theory/covariances_fourier_nobin_simple.c"

#include "init_LSSxCMB.c"

#include "../cosmolike_core/cfftlog/utils_complex.h"
// Naming convention:
// l = galaxy positions ("l" as in "lens sample")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "source sample")
// And alphabetical order

// used in checking existence of output directory
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//

void run_cov_AB_CD(char ABCD[2][4], char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);

void run_cov_AB_CD(char ABCD[2][4], char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2, int start){
  int nl1,nl2,i,j,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  int n12[2], z_ar[4], is_ls[4], N_start[2];
  n12[0] = n1; n12[1] = n2;
  for(i=0;i<4;i++) {
    is_ls[i]=0; //l->1, s-> -1, others ->0 (default)
  }

  for(i=0;i<2;i++){
    if(strcmp(ABCD[i],"ss")==0){
      z_ar[2*i] = Z1(n12[i]); z_ar[2*i+1] = Z2(n12[i]);
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i]=-1; is_ls[2*i+1]=-1;
      N_start[i] = like.Ncl*n12[i];
    }else if(strcmp(ABCD[i],"ls")==0){
      z_ar[2*i] = ZL(n12[i]); z_ar[2*i+1] = ZS(n12[i]);
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i]=1; is_ls[2*i+1]=-1; // l->1, s-> -1, others ->0 (default)
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+n12[i]);
    }else if(strcmp(ABCD[i],"ll")==0){
      z_ar[2*i]=n12[i]; z_ar[2*i+1]=n12[i];
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i]=1; is_ls[2*i+1]=1;
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n12[i]);
    }else if(strcmp(ABCD[i],"lk")==0){
      z_ar[2*i] = n12[i]; z_ar[2*i+1] = -1; // z_ar = -1 for cmb kappa field (for no z bin index)
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i]=1;
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n12[i]);
    }else if(strcmp(ABCD[i],"ks")==0){
      z_ar[2*i] = -1; z_ar[2*i+1] = n12[i];
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i+1]=-1;
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n12[i]);
    }else if(strcmp(ABCD[i],"kk")==0){
      z_ar[2*i] = -1; z_ar[2*i+1] = -1;
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin);
    }else if(strcmp(ABCD[i],"ly")==0){
      z_ar[2*i] = n12[i]; z_ar[2*i+1] = -2; // z_ar = -2 for y
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i]=1;
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin+1+n12[i]);
    }else if(strcmp(ABCD[i],"sy")==0){
      z_ar[2*i] = n12[i]; z_ar[2*i+1] = -2;
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      is_ls[2*i]=-1;
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin+1+tomo.clustering_Nbin+n12[i]);
    }else if(strcmp(ABCD[i],"ky")==0){
      z_ar[2*i] = -1; z_ar[2*i+1] = -2;
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin+1+tomo.clustering_Nbin+tomo.shear_Nbin);
    }else if(strcmp(ABCD[i],"yy")==0){
      z_ar[2*i] = -2; z_ar[2*i+1] = -2;
      printf("ABCD[%d]=%s, bins = %d (%d, %d)\n", i,ABCD[i], n12[i],z_ar[2*i],z_ar[2*i+1]);
      N_start[i] = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin+1+tomo.clustering_Nbin+tomo.shear_Nbin+1);
    }
  }

  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = 1.;
      for(i=0; i<2; i++) {
        if(is_ls[i]==1) {weight *= test_kmax(ell[nl1],z_ar[i]);}
      }
      for(i=2; i<4; i++) {
        if(is_ls[i]==1) {weight *= test_kmax(ell[nl2],z_ar[i]);}
      }

      if(weight){
        if (nl1 == nl2){c_g = cov_G_AB_CD(ABCD, ell[nl1],dell[nl1],z_ar, is_ls);}
        if (covparams.ng){c_ng = cov_NG_AB_CD(ABCD, ell[nl1],ell[nl2],z_ar, is_ls);}
      }

      i=N_start[0]+nl1;
      j=N_start[1]+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],z_ar[0],z_ar[1],z_ar[2],z_ar[3],c_g,c_ng);
    }
  }
  fclose(F1);
}



#define N_scenarios 7
int main(int argc, char** argv)
{
  
  int i,l,m,n,o,s,p,nl1,t,k;
  char OUTFILE[400],filename[400],arg1[400],arg2[400];
  

  // Y1 corresponds to DESC SRD Y1, Y6 corresponds to assuming that we cover the full SO area=0.4*fsky and at a depth of 26.1 which is in a range of reasonable scenarios (see https://github.com/LSSTDESC/ObsStrat/tree/static/static )
  // Roman fiducial: 2000deg^2, nsrc=51, nlens=66
  // Roman Wide: 18000deg^2, nsrc=43, nlens=50
  // Roman Wide2: 10000deg^2, nsrc=43, nlens=50
  double area_table[N_scenarios]={12300.0,16500.0,18000.,2000.,2000.,18000.,10000.}; 
  // double nsource_table[3]={11.0,23.0,28.0};
  // double nlens_table[3]={18.0,41.0,48.0};

  double nsource_table[N_scenarios]={10.7,22.5,27.1, 51.0,51.0, 43.0,43.0};
  double nlens_table[N_scenarios]={13.1,26.8,32.0, 66.0,66.0, 50.0,50.0};

  // Roman uses optimistic scenario (Table 2 of 2004.05271)
  double sigma_zphot_shear[N_scenarios]={0.05,0.05,0.05, 0.01,0.01, 0.02,0.02};
  double sigma_zphot_clustering[N_scenarios]={0.03,0.03,0.03, 0.01,0.01, 0.02,0.02};
  
  char survey_designation[N_scenarios][200]={"LSSTxSO_Y1","LSSTxSO_Y6","LSSTxSO_Y10", "RomanxSO","RomanxS4","RomanWidexS4","RomanWide2xS4"};
  
  char source_zfile[N_scenarios][400]={"src_LSSTY1","src_LSSTY6","src_LSSTY10", \
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm"};

#ifdef ONESAMPLE
  char lens_zfile[N_scenarios][400]={"src_LSSTY1","src_LSSTY6","src_LSSTY10", \
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_lensing_fine_bin_norm"};
  for (i=0; i<N_scenarios; i++) {nlens_table[i] = nsource_table[i];}

#else
  char lens_zfile[N_scenarios][400]={"lens_LSSTY1","lens_LSSTY6","lens_LSSTY10", \
                  "zdistri_WFIRST_LSST_clustering_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_clustering_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_clustering_fine_bin_norm",\
                  "zdistri_WFIRST_LSST_clustering_fine_bin_norm"};
#endif

  int hit=atoi(argv[1]);
  Ntable.N_a=100;
  k=1;
  
  t = atoi(argv[2]);
  
  //RUN MODE setup
  // init_cosmo_runmode("halofit");
  init_cosmo_runmode("halomodel");
  // init_binning_fourier(20,30.0,3000.0,3000.0,21.0,10,10);
  // init_binning_fourier(15,20.0,3000.0,3000.0,0.0,10,10);
  init_binning_fourier(25,20.0,7979.0,7979.0,0.0,10,10);
  init_survey(survey_designation[t],nsource_table[t],nlens_table[t],area_table[t]);
  // set sigma_zphot_shear and sigma_zphot_clustering for gaussian nz error
  init_priors(0.002,sigma_zphot_shear[t],0.001,0.001,sigma_zphot_clustering[t],0.001,0.001);

  sprintf(arg1,"zdistris/%s",source_zfile[t]);
  sprintf(arg2,"zdistris/%s",lens_zfile[t]); 
#ifdef ONESAMPLE
  init_galaxies(arg1,arg2,"gaussian","gaussian","source_std","lens=src");
#else
  if(t < 3) { // LSST
    init_galaxies(arg1,arg2,"gaussian","gaussian","source_std","LSST_gold");
  } else { // Roman
    init_galaxies(arg1,arg2,"gaussian","gaussian","source_std","WF_SN10");
  }
#endif
  init_IA("none","GAMA");
  char probe[] = "10x2pt";
  init_probes(probe);

  if(t==0) init_cmb("so_Y1");
  if( (t>0) || (t<=3) ) init_cmb("so_Y5");
  if(t==4 || t==5 || t==6) init_cmb("s4");

  //set l-bins for shear, ggl, clustering, clusterWL
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  double *ell, *dell, *ell_Cluster, *dell_Cluster;
  ell=create_double_vector(0,like.Ncl-1);
  dell=create_double_vector(0,like.Ncl-1);
  int j=0;
  for(i=0;i<like.Ncl;i++){
    ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
    dell[i]=exp(log(like.lmin)+(i+1)*logdl)-exp(log(like.lmin)+(i*logdl));
    if(ell[i]<like.lmax_shear) printf("%le\n",ell[i]);
  } 

  covparams.ng = 1;
  covparams.cng = 1;
#ifdef SSCONLY
  covparams.cng = 0;
#endif

  printf("----------------------------------\n");  
  sprintf(survey.name,"%s_area%le_ng%le_nl%le",survey_designation[t],survey.area,survey.n_gal,survey.n_lens);
  printf("area: %le n_source: %le n_lens: %le\n",survey.area,survey.n_gal,survey.n_lens);

  // Setup output dir
  sprintf(covparams.outdir,"out_cov_%s_%s", survey_designation[t], probe);
#ifdef ONESAMPLE
  sprintf(covparams.outdir,"%s_1sample", covparams.outdir);
#endif
#ifdef SLOW
  sprintf(covparams.outdir,"%s_slow", covparams.outdir);
#endif
#ifdef RUN_FFT
  sprintf(covparams.outdir,"%s_fft", covparams.outdir);
#endif
#ifdef SSCONLY
  sprintf(covparams.outdir,"%s_ssc", covparams.outdir);
#endif

  sprintf(covparams.outdir,"%s/", covparams.outdir);

  // if outdir doesn't exist, create one
  struct stat st = {0};
  if (stat(covparams.outdir, &st) == -1) {
      mkdir(covparams.outdir, 0700);
  }

  printf("----------------------------------\n");
  char ABCD[2][4];
  const char *probes[10] = {"ss", "ls", "ll", "lk", "ks", "kk",\
                            "ly", "sy", "ky", "yy"};
  const int Npower[10] = {tomo.shear_Npowerspectra,tomo.ggl_Npowerspectra,tomo.clustering_Npowerspectra,\
                          tomo.clustering_Nbin,tomo.shear_Nbin,1,\
                          tomo.clustering_Nbin,tomo.shear_Nbin,1,1};
  int N_blocks[56], i_block=1; // There are 55 blocks in total
  int i_probe1[56], i_probe2[56];
  N_blocks[0]=0;
  printf("Block index starts from 1:\n");
  for(i=0;i<10;i++){
    for(j=0;j<=i;j++){
      if(i==j){
        N_blocks[i_block] = N_blocks[i_block-1] + Npower[i]*(Npower[i]+1)/2;
      }else{
        N_blocks[i_block] = N_blocks[i_block-1] + Npower[i]*Npower[j];
      }
      printf("%s-%s ends at %d\n", probes[i], probes[j], N_blocks[i_block]);
      i_probe1[i_block] = i;
      i_probe2[i_block] = j;
      i_block++;
    }
  }

  if(hit>N_blocks[55]){
    printf("Cov index %d exceeding the total number of blocks N=%d\n", hit, N_blocks[55]);
    exit(1);
  }

  i_block=1;
  while(hit>N_blocks[i_block]){ // until hit<=N_blocks[i_block]
    i_block++;
  }
  int m0, i_pro1, i_pro2;
  i_pro1 = i_probe1[i_block]; i_pro2 = i_probe2[i_block];
  sprintf(OUTFILE,"%s_%s%s_cov_Ncl%d_Ntomo%d",survey.name,probes[i_pro1],probes[i_pro2], like.Ncl,tomo.shear_Nbin);
  hit -= N_blocks[i_block-1];

  for (l=0;l<Npower[i_pro1]; l++){
    m0 = (i_pro1==i_pro2) ? l : 0;
    for (m=m0;m<Npower[i_pro2]; m++){
      if(k==hit){
        sprintf(ABCD[0],"%s",probes[i_pro1]);
        sprintf(ABCD[1],"%s",probes[i_pro2]);
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_AB_CD(ABCD,OUTFILE,covparams.outdir,ell,dell,l,m,k+N_blocks[i_block-1]);
      }
      k=k+1;
    }
  }

  printf("number of cov blocks for parallelization: %d\n",N_blocks[55]); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

