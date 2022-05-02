#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

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
#include "../cosmolike_core/theory/CMBxLSS_fourier.c"
// #include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/init_baryon.c"
#include "init_LSSxCMB.c"

#include "../cosmolike_core/cfftlog/utils_complex.h"
// Naming convention:
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")
// And alphabetical order

typedef double (*C_tomo_pointer)(double l, int n1, int n2);
void twopoint_via_hankel(double **xi, double *logthetamin, double *logthetamax, C_tomo_pointer C_tomo, int ni, int nj, int N_Bessel);



typedef struct input_cosmo_params_local {
    double omega_m;
    double sigma_8;
    double n_s;
    double w0;
    double wa;
    double omega_b;
    double h0;
    double MGSigma;
    double MGmu;
} input_cosmo_params_local;

typedef struct input_nuisance_params_local {
    double bias[10];
    double source_z_bias[10];
    double source_z_s;
    double lens_z_bias[10];
    double lens_z_s;
    double shear_m[10];
    double A_ia;
    double eta_ia;
#ifdef TEST_CALIB
    double gas[15];
#else
    double gas[11];
#endif
} input_nuisance_params_local;

void print_cosmo_params(input_cosmo_params_local ic);
void print_nuisance_params(input_nuisance_params_local in);

double C_shear_tomo_sys(double ell,int z1,int z2);
// double C_cgl_tomo_sys(double ell_Cluster,int zl,int nN, int zs);
double C_gl_tomo_sys(double ell,int zl,int zs);
double C_ks_sys(double ell, int zs);
void set_data_shear(double *ell, double *data, int start);
void set_data_ggl(double *ell, double *data, int start);
void set_data_clustering(double *ell, double *data, int start);
void set_data_gk(double *ell, double *data, int start);
void set_data_ks(double *ell, double *data, int start);
void set_data_kk(double *ell, double *data, int start);
void set_data_gy(double *ell, double *data, int start);
void set_data_sy(double *ell, double *data, int start);
void set_data_ky(double *ell, double *data, int start);
void set_data_yy(double *ell, double *data, int start);
// void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q,double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope, double Q1, double Q2, double Q3);
// double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q,double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope, double Q1, double Q2, double Q3);
void compute_data_vector(char *details, input_cosmo_params_local ic, input_nuisance_params_local in);

double log_multi_like(input_cosmo_params_local ic, input_nuisance_params_local in);

int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);
int get_N_ggl(void);
int get_N_ell(void);


int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}
int get_N_ggl(void){
  return tomo.ggl_Npowerspectra;
}
int get_N_ell(void){
  return like.Ncl;
}

double C_shear_tomo_sys(double ell, int z1, int z2)
{
  double C;
  // C= C_shear_tomo_nointerp(ell,z1,z2);
  // if(like.IA==1) C+=C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  
  // if(like.IA!=1) C= C_shear_tomo_nointerp(ell,z1,z2);
  if(like.IA==4) C= C_shear_shear_IA(ell,z1,z2);
  else{printf("IA not supported in this analysis.\n"); exit(1);}
  // clock_t t1, t2; t1=clock();
  // if(like.IA==1) C = C_shear_tomo_nointerp(ell,z1,z2)+C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  // t2 = clock(); printf("shear %d-%d: %le\n",z1,z2, (double)(t2-t1)/CLOCKS_PER_SEC);
  // if(like.IA==2) C += C_II_lin_nointerp(ell,z1,z2)+C_GI_lin_nointerp(ell,z1,z2);  
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
  //printf("%le %d %d %le\n",ell,z1,z2,C_shear_tomo_nointerp(ell,z1,z2)+C_II_JB_nointerp(ell,z1,z2)+C_GI_JB_nointerp(ell,z1,z2));
return C;
}

double C_gl_tomo_sys(double ell,int zl,int zs)
{
  double C;
  // C=C_gl_tomo_nointerp(ell,zl,zs); 
  // if(like.IA==1) C += C_gI_nointerp(ell,zl,zs);
  
  // if(like.IA!=1) C=C_gl_tomo_nointerp(ell,zl,zs);
  if(like.IA==4) C = C_ggl_IA(ell,zl,zs);
  else{printf("IA not supported in this analysis.\n"); exit(1);}
  // if(like.IA==2) C += C_gI_lin_nointerp(ell,zl,zs);
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
return C;
}

double C_ks_sys(double ell, int zs)
{
   double C;
   C = C_ks_nointerp(ell,zs);
   if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
   return C;
}

double C_sy_sys(double ell, int zs)
{
   double C;
   C = C_sy_nointerp(ell,zs);
   if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
   return C;
}

void set_data_shear(double *ell, double *data, int start)
{
  int i,z1,z2,nz;
  double a;
  int Ncl = like.Ncl;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < Ncl; i++){
      if (ell[i] < like.lmax_shear){ data[Ncl*nz+i] = C_shear_tomo_sys(ell[i],z1,z2);}
      else {data[Ncl*nz+i] = 0.;}
    }
  }
}

void set_data_ggl(double *ell, double *data, int start)
{
  int i, zl,zs,nz; 
  int Ncl = like.Ncl;
  for (nz = 0; nz < tomo.ggl_Npowerspectra; nz++){
    zl = ZL(nz); zs = ZS(nz);
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],zl) && ell[i] < like.lmax_shear){
        data[start+(Ncl*nz)+i] = C_gl_tomo_sys(ell[i],zl,zs);
      }
      else{
        data[start+(Ncl*nz)+i] = 0.;
      }
    }
  }
}

void set_data_clustering(double *ell, double *data, int start){
  int i, nz;
  int Ncl = like.Ncl;
  for (nz = 0; nz < tomo.clustering_Npowerspectra; nz++){
    //printf("%d %e %e\n",nz, gbias.b[nz][1],pf_photoz(gbias.b[nz][1],nz));
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],nz)){data[start+(Ncl*nz)+i] = C_cl_tomo_nointerp(ell[i],nz,nz);}
      else{data[start+(Ncl*nz)+i] = 0.;}
      //printf("%d %d %le %le\n",nz,nz,ell[i],data[Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + nz)+i]);
    }
  }
}


void set_data_gk(double *ell, double *data, int start)
{
  for (int nz=0; nz<tomo.clustering_Nbin; nz++){
    for (int i=0; i<like.Ncl; i++){
       if (ell[i]<like.lmax_kappacmb && ell[i]>like.lmin_kappacmb && test_kmax(ell[i],nz)){
          data[start+(like.Ncl*nz)+i] = C_gk_nointerp(ell[i],nz);
       }
       else{
          data[start+(like.Ncl*nz)+i] = 0.;
       }
    }
  }
}

void set_data_ks(double *ell, double *data, int start)
{
   for (int nz=0; nz<tomo.shear_Nbin; nz++){
      for (int i=0; i<like.Ncl; i++){
         if (ell[i]<like.lmax_kappacmb && ell[i]>like.lmin_kappacmb && ell[i] < like.lmax_shear) {
            data[start+(like.Ncl*nz)+i] = C_ks_sys(ell[i],nz);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      }
   }
}

void set_data_kk(double *ell, double *data, int start)
{
   for (int i=0; i<like.Ncl; i++){
      if (ell[i]<like.lmax_kappacmb && ell[i]>like.lmin_kappacmb){
         data[start+i] = C_kk_nointerp(ell[i]);
      }
      else{
         data[start+i] = 0.;
      }
   }
}

void set_data_gy(double *ell, double *data, int start)
{
  for (int nz=0; nz<tomo.clustering_Nbin; nz++){
    for (int i=0; i<like.Ncl; i++){
       if (ell[i]<like.lmax_y && ell[i]>like.lmin_y && test_kmax(ell[i],nz)){
          data[start+(like.Ncl*nz)+i] = C_gy_nointerp(ell[i],nz);
       }
       else{
          data[start+(like.Ncl*nz)+i] = 0.;
       }
    }
  }
}

void set_data_sy(double *ell, double *data, int start)
{
   for (int nz=0; nz<tomo.shear_Nbin; nz++){
      for (int i=0; i<like.Ncl; i++){
         if (ell[i]<like.lmax_y && ell[i]>like.lmin_y && ell[i] < like.lmax_shear) {
            data[start+(like.Ncl*nz)+i] = C_sy_sys(ell[i],nz);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      }
   }
}

void set_data_ky(double *ell, double *data, int start)
{
   for (int i=0; i<like.Ncl; i++){
      if (ell[i]<like.lmax_y && ell[i]>like.lmin_y && ell[i]<like.lmax_kappacmb && ell[i]>like.lmin_kappacmb){
         data[start+i] = C_ky_nointerp(ell[i]);
      }
      else{
         data[start+i] = 0.;
      }
   }
}

void set_data_yy(double *ell, double *data, int start)
{
   for (int i=0; i<like.Ncl; i++){
      if (ell[i]<like.lmax_y && ell[i]>like.lmin_y){
         data[start+i] = C_yy_nointerp(ell[i]);
      }
      else{
         data[start+i] = 0.;
      }
   }
}

int set_cosmology_params(input_cosmo_params_local ic)
{
  cosmology.Omega_m=ic.omega_m;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  cosmology.sigma_8=ic.sigma_8;
  cosmology.n_spec= ic.n_s;
  cosmology.w0=ic.w0;
  cosmology.wa=ic.wa;
  cosmology.omb=ic.omega_b;
  cosmology.h0=ic.h0;
  cosmology.MGSigma=ic.MGSigma;
  cosmology.MGmu=ic.MGmu;

  if (cosmology.Omega_m < 0.1 || cosmology.Omega_m > 0.6) return 0;
  if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
  if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
  if (cosmology.n_spec < 0.85 || cosmology.n_spec > 1.05) return 0;
  if (cosmology.w0 < -2.0 || cosmology.w0 > -0.0) return 0;
  if (cosmology.wa < -2.5 || cosmology.wa > 2.5) return 0;
  if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  if (cosmology.MGmu < -3.0 || cosmology.MGmu > 3.0) return 0;
  if (cosmology.MGSigma < -3.0 || cosmology.MGSigma > 3.0) return 0; // DESY1 extension paper flat priors

  return 1;
}

void set_nuisance_shear_calib(double *M)
{
  for(int i=0;i<tomo.shear_Nbin;i++) {nuisance.shear_calibration_m[i] = M[i];}
}

int set_nuisance_shear_photoz(double *SP, double SPS1)
{
  for(int i=0;i<tomo.shear_Nbin;i++) {
    nuisance.bias_zphot_shear[i] = SP[i];
    nuisance.sigma_zphot_shear[i]=SPS1;
    if (nuisance.sigma_zphot_shear[i]<0.0001) return 0;
  }
  return 1;
}

int set_nuisance_clustering_photoz(double *CP, double CPS1)
{
  for(int i=0;i<tomo.clustering_Nbin;i++) {
    nuisance.bias_zphot_clustering[i]=CP[i];
    nuisance.sigma_zphot_clustering[i]=CPS1;
    if (nuisance.sigma_zphot_clustering[i]<0.0001) return 0;
  }
  return 1;
}

int set_nuisance_ia(double A_ia, double eta_ia)
{
  nuisance.A_ia=A_ia;  
  nuisance.eta_ia=eta_ia;
  if (nuisance.A_ia < -5.0 || nuisance.A_ia > 5.0) return 0;
  if (nuisance.eta_ia < -5.0 || nuisance.eta_ia> 5.0) return 0;
  return 1;
}

int set_nuisance_gbias(double *B)
{
  for(int i=0;i<tomo.clustering_Nbin;i++) {
    gbias.b[i] = B[i];
#ifdef ONESAMPLE
    if (gbias.b[i] < 0.4 || gbias.b[i] > 5.0) return 0;
#else
    if (gbias.b[i] < 0.4 || gbias.b[i] > 3.0) return 0;
#endif
  }
  return 1;
} 

int set_nuisance_gas(double *p_gas)
{
   nuisance.gas_Gamma_KS=p_gas[0];
   nuisance.gas_beta=p_gas[1];
   nuisance.gas_lgM0=p_gas[2];
   nuisance.gas_alpha=p_gas[3];
   nuisance.gas_A_star=p_gas[4];
   nuisance.gas_lgM_star=p_gas[5];
   nuisance.gas_sigma_star=p_gas[6];
   nuisance.gas_lgT_w=p_gas[7];
   nuisance.gas_f_H=p_gas[8];
   nuisance.gas_eps1=p_gas[9];
   nuisance.gas_eps2=p_gas[10];
   if (nuisance.gas_Gamma_KS < 1.05 || nuisance.gas_Gamma_KS > 1.35) return 0;
   if (nuisance.gas_beta < 0.2 || nuisance.gas_beta > 1.0) return 0;
   if (nuisance.gas_lgM0 < 12.5 || nuisance.gas_lgM0 > 15.0) return 0;
   if (nuisance.gas_alpha < 0.5 || nuisance.gas_alpha > 1.5) return 0;

   if (nuisance.gas_lgT_w < 6.0 || nuisance.gas_lgT_w > 7.0) return 0;
   if (nuisance.gas_eps1 < -0.3 || nuisance.gas_eps1 > 0.3) return 0;
   if (nuisance.gas_eps2 < -0.3 || nuisance.gas_eps2 > 0.3) return 0;
   if (nuisance.gas_f_H < 0.7 || nuisance.gas_f_H > 0.8) return 0;

   // need to check the ranges below if they are varied
   if (nuisance.gas_lgM_star < 12.0 || nuisance.gas_lgM_star > 14.0) return 0;
   if (nuisance.gas_sigma_star < 1.0 || nuisance.gas_sigma_star > 1.5) return 0;
   if (nuisance.gas_A_star < 0.02 || nuisance.gas_A_star > 0.04) return 0;

#ifdef TEST_CALIB
   nuisance.gas_beta_v2=p_gas[11];
   nuisance.gas_lgM0_v2=p_gas[12];
   nuisance.gas_eps1_v2=p_gas[13];
   nuisance.gas_eps2_v2=p_gas[14];

   if (nuisance.gas_beta_v2 < 0.2 || nuisance.gas_beta_v2 > 1.0) return 0;
   if (nuisance.gas_lgM0_v2 < 12.5 || nuisance.gas_lgM0_v2 > 15.0) return 0;
   if (nuisance.gas_eps1_v2 < -0.3 || nuisance.gas_eps1_v2 > 0.3) return 0;
   if (nuisance.gas_eps2_v2 < -0.3 || nuisance.gas_eps2_v2 > 0.3) return 0;

#endif


return 1;
}

// double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu,\
//                        double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10,\
//                        double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1,\
//                        double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1,\
//                        double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10,\
//                        double A_ia, double eta_ia,\
//                        double Gamma_KS, double beta, double lgM0, double alpha,\
//                        double A_star, double lgM_star, double sigma_star,\
//                        double lgT_w, double f_H, double eps1, double eps2)
double log_multi_like(input_cosmo_params_local ic, input_nuisance_params_local in)
{
  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0, log_L=0.0;;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
  }

  if (set_cosmology_params(ic)==0){
    printf("Cosmology out of bounds\n");
    return -1.0e15;
  }
  set_nuisance_shear_calib(in.shear_m);
  if (set_nuisance_shear_photoz(in.source_z_bias, in.source_z_s)==0){
    printf("Shear photo-z sigma too small\n");
    return -1.0e15;
  }
  if (set_nuisance_clustering_photoz(in.lens_z_bias, in.lens_z_s)==0){
    printf("Clustering photo-z sigma too small\n");
    return -1.0e15;
  }
  if (set_nuisance_ia(in.A_ia,in.eta_ia)==0){
    printf("IA parameters out of bounds\n");
    return -1.0e15; 
  }
  if (set_nuisance_gbias(in.bias)==0){
    printf("Bias out of bounds\n");
    return -1.0e15;
  }
  if (set_nuisance_gas(in.gas)==0){
    printf("Gas parameters out of bounds\n");
    return -1.0e15;
  }
       
  printf("like: %le %le %le %le %le %le %le %le\n",cosmology.Omega_m, cosmology.Omega_v,cosmology.sigma_8,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,cosmology.h0); 
  // printf("like %le %le %le %le\n",gbias.b[0][0], gbias.b[1][0], gbias.b[2][0], gbias.b[3][0]);    
  // for (i=0; i<10; i++){
  //   printf("nuisance %le %le %le\n",nuisance.shear_calibration_m[i],nuisance.bias_zphot_shear[i],nuisance.sigma_zphot_shear[i]);
  // }

  log_L_prior=0.0;
  // if(like.Aubourg_Planck_BAO_SN==1) log_L_prior+=log_L_Planck_BAO_SN();
  // if(like.SN==1) log_L_prior+=log_L_SN();
  //if(like.BAO==1) log_L_prior+=log_L_BAO();
  // if(like.Planck==1) log_L_prior+=log_L_Planck();
  // if(like.Planck15_BAO_w0wa==1) log_L_prior+=log_L_Planck15_BAO_w0wa();//CH
  //if(like.Planck15_BAO_H070p6_JLA_w0wa==1) log_L_prior+=log_L_Planck15_BAO_H070p6_JLA_w0wa();//CH
  // if(like.IA!=0) log_L_prior+=log_L_ia();
  // if(like.IA!=0) log_L_prior+=log_like_f_red();
  if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  if(like.clphotoz!=0) log_L_prior+=log_L_clphotoz();
  if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();
  // if(like.IA!=0) {
  //   log_L = 0.0;
  //   log_L -= pow((nuisance.A_ia - prior.A_ia[0])/prior.A_ia[1],2.0);
  //   log_L -= pow((nuisance.eta_ia - prior.eta_ia[0])/prior.eta_ia[1],2.0);
  //   log_L_prior+=0.5*log_L;
  // }
 
  // if(like.clusterMobs==1) log_L_prior+=log_L_clusterMobs();
 
  // printf("%d %d %d %d\n",like.BAO,like.wlphotoz,like.clphotoz,like.shearcalib);
  // printf("logl %le %le %le %le\n",log_L_shear_calib(),log_L_wlphotoz(),log_L_clphotoz(),log_L_clusterMobs());
  int start=0;  
  // clock_t t1, t2;
  // t1 = clock();
  if(like.shear_shear==1) {
    set_data_shear(ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  // t2 = clock(); printf("shear: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if(like.shear_pos==1){
    set_data_ggl(ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  // t2 = clock(); printf("ggl: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if(like.pos_pos==1){
    set_data_clustering(ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }
  // t2 = clock(); printf("gg: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if(like.gk==1) {
    set_data_gk(ell, pred, start);
    start += like.Ncl*tomo.clustering_Nbin;
  }
  // t2 = clock(); printf("gk: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if(like.ks==1) {
    set_data_ks(ell, pred, start);
    start += like.Ncl*tomo.shear_Nbin;
  }
  // t2 = clock(); printf("ks: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if(like.kk==1) {
    set_data_kk(ell, pred, start);
    start += like.Ncl;
  }
  // t2 = clock(); printf("kk: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if (like.gy) {
    set_data_gy(ell, pred, start);
    start += like.Ncl * tomo.clustering_Nbin;
  }
  if (like.sy) {
    set_data_sy(ell, pred, start);
    start += like.Ncl * tomo.shear_Nbin;
  }
  if (like.ky) {
    set_data_ky(ell, pred, start);
    start += like.Ncl;
  }
  if (like.yy) {
    set_data_yy(ell, pred, start);
    start += like.Ncl;
  }

  chisqr=0.0;
  for (i=0; i<like.Ndata; i++){
    if(dvmask_read(1,i)){
      for (j=0; j<like.Ndata; j++){
        if(dvmask_read(1,j)){
          a=(pred[i]-data_read(1,i))*invcov_read(1,i,j)*(pred[j]-data_read(1,j));
          chisqr += a;
        }
      }
    }
    // if (fabs(data_read(1,i)) < 1.e-25){
    //    printf("%d %le %le %le\n",i,data_read(1,i),pred[i],invcov_read(1,i,i));
    // }
  }
  // t2 = clock(); printf("chisqr: %le\n", (double)(t2-t1)/CLOCKS_PER_SEC); t1 = t2;
  if (chisqr<0.0){
    printf("error: chisqr = %le\n",chisqr);
    //exit(EXIT_FAILURE);
  }
  if (isnan(chisqr)){
    print_cosmo_params(ic);
    print_nuisance_params(in);
  }
  // printf("chisqr, log_L_prior: %le, %le\n",chisqr, log_L_prior);
  log_L = -0.5*chisqr+log_L_prior;
  printf("log_L: %le\n",log_L);
  if(log_L < -1.0e15 || isnan(log_L)) {return -1.0e15;}
  return log_L;
}

void print_cosmo_params(input_cosmo_params_local ic){
  double p_cosmo[9]={ic.omega_m,ic.sigma_8,ic.n_s,ic.w0,ic.wa,ic.omega_b,ic.h0,ic.MGSigma,ic.MGmu};
  for(int i=0; i<9; i++){ printf("%le, ", p_cosmo[i]); }
  printf("\n");
}

void print_nuisance_params(input_nuisance_params_local in){
  printf("gbias: ");
  for(int i=0; i<tomo.clustering_Nbin; i++){ printf("%le, ", in.bias[i]); }
  printf("\n"); printf("dz_src: ");
  for(int i=0; i<tomo.shear_Nbin; i++){ printf("%le, ", in.source_z_bias[i]); }
  printf("(sig_zsrc) %le, ", in.source_z_s); 
  printf("\n"); printf("dz_lens: ");
  for(int i=0; i<tomo.clustering_Nbin; i++){ printf("%le, ", in.lens_z_bias[i]); }
  printf("(sig_zlens) %le, ", in.lens_z_s);
  printf("\n"); printf("m: ");
  for(int i=0; i<tomo.shear_Nbin; i++){ printf("%le, ", in.source_z_bias[i]); }
  printf("\n");
  printf("IA (A, eta): %le, %le\n", in.A_ia, in.eta_ia);
  printf("gas: ");
#ifdef TEST_CALIB
  int N_gaspara = 15;
#else
  int N_gaspara = 11;
#endif
  for(int i=0; i<N_gaspara; i++){ printf("%le, ", in.gas[i]); }
  printf("\n");
}

void compute_data_vector(char *details, input_cosmo_params_local ic, input_nuisance_params_local in)
{

  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  // static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
  }
// for (l=0;l<like.Ncl;l++){
//   printf("%d %le\n",i,ell[l]);
// }
  // clock_t t1, t2;

  // t1 = clock();
  
  set_cosmology_params(ic);
  set_nuisance_shear_calib(in.shear_m);
  set_nuisance_shear_photoz(in.source_z_bias, in.source_z_s);
  set_nuisance_clustering_photoz(in.lens_z_bias, in.lens_z_s);
  set_nuisance_ia(in.A_ia,in.eta_ia);
  set_nuisance_gbias(in.bias);
  set_nuisance_gas(in.gas);

  int start=0;  
  if(like.shear_shear==1) {
    set_data_shear(ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  if(like.shear_pos==1){
    printf("ggl, start %d\n", start);
    set_data_ggl(ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  }
  if(like.pos_pos==1){
    printf("clustering, start %d\n", start);
    set_data_clustering(ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }

  if(like.gk==1) {
    printf("Computing data vector: gk, start %d\n", start);
    set_data_gk(ell, pred, start);
    start += like.Ncl * tomo.clustering_Nbin;
  }
  if(like.ks==1) {
    printf("Computing data vector: ks, start %d\n", start);
    set_data_ks(ell, pred, start);
    start += like.Ncl * tomo.shear_Nbin;
  }
  if (like.kk) {
    printf("Computing data vector: kk, start %d\n", start);
    set_data_kk(ell, pred, start);
    start += like.Ncl;
  }

  if (like.gy) {
    printf("Computing data vector: gy, start %d\n", start);
    set_data_gy(ell, pred, start);
    start += like.Ncl * tomo.clustering_Nbin;
  }
  if (like.sy) {
    printf("Computing data vector: sy, start %d\n", start);
    set_data_sy(ell, pred, start);
    start += like.Ncl * tomo.shear_Nbin;
  }
  if (like.ky) {
    printf("Computing data vector: ky, start %d\n", start);
    set_data_ky(ell, pred, start);
    start += like.Ncl;
  }
  if (like.yy) {
    printf("Computing data vector: yy, start %d\n", start);
    set_data_yy(ell, pred, start);
    start += like.Ncl;
  }

  FILE *F;
  char filename[300];
  if (strstr(details,"FM") != NULL){
    sprintf(filename,"%s",details);
  }
  else {sprintf(filename,"datav/%s_%s",like.probes,details);}
  F=fopen(filename,"w");
  for (i=0;i<like.Ndata; i++){  
    fprintf(F,"%d %le\n",i,pred[i]);
    //printf("%d %le\n",i,pred[i]);
  }
  fclose(F);

}

void save_zdistr_sources(int zs){
  double z,dz =(redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) for source redshift bin %d\n",zs);
  
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_sources_bin%d.txt",zs);
   F1 = fopen(filename,"w");
   for (z =redshift.shear_zdistrpar_zmin; z< redshift.shear_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, zdistr_photoz(z,zs));
   }
}


void save_zdistr_lenses(int zl){
   double z,dz =(redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) and bias b(z) for lens redshift bin %d\n",zl);
   
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_lenses_bin%d.txt", zl);
   F1 = fopen(filename,"w");
   for (z =redshift.clustering_zdistrpar_zmin; z< redshift.clustering_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, pf_photoz(z,zl));
   }
}


 int main(int argc, char** argv)
{
  int i;
  char arg1[400],arg2[400],arg3[400];
/* here, do your time-consuming job */
  int sce=atoi(argv[1]);

  int N_scenarios=2;

  double sigma_zphot_shear[3]={0.05,0.05};
  double sigma_zphot_clustering[3]={0.03,0.03};

  double area_table[2]={12300.0,16500.0}; // Y1 corresponds to DESC SRD Y1, Y6 corresponds to assuming that we cover the full SO area=0.4*fsky and at a depth of 26.1 which is in a range of reasonable scenarios (see https://github.com/LSSTDESC/ObsStrat/tree/static/static )
  // double nsource_table[2]={11.0,23.0};
  // double nlens_table[2]={18.0,41.0};

  double nsource_table[2]={10.7,22.5};
  double nlens_table[2]={13.1,26.8};

  char survey_designation[2][200]={"LSSTxSO_Y1","LSSTxSO_Y6"};
  char tomo_binning_source[2][200]={"source_std","source_std"};
  // even for lens=src, this setting is valid, because the input lens/src zfiles are generated with the same z-cut:(0.2,1.2)

  char source_zfile[2][400]={"src_LSSTY1","src_LSSTY6"};
#ifdef ONESAMPLE
  char lens_zfile[2][400]={"src_LSSTY1","src_LSSTY6"};
  char tomo_binning_lens[2][200]={"lens=src","lens=src"};
  nlens_table[0] = nsource_table[0];
  nlens_table[1] = nsource_table[1];
  sigma_zphot_clustering[0] = sigma_zphot_shear[0];
  sigma_zphot_clustering[1] = sigma_zphot_shear[1];
#else
  char lens_zfile[2][400]={"lens_LSSTY1", "lens_LSSTY6"};
  char tomo_binning_lens[2][200]={"LSST_gold","LSST_gold"};
#endif
  char cmb_yr[2][100]={"so_Y1","so_Y5"};

  init_cosmo_runmode("halomodel");
  // init_cosmo_runmode("halofit");
  // double z_max_limit = 4.;
  // limits.a_min = 1./(1.+z_max_limit);
  // limits.a_min_hm = 1./(1.+10);

  // init_bary(argv[2]);
  // init_binning_fourier(15,20.0,3000.0,3000.0,21.0,10,10);
  init_binning_fourier(25,20.0,7979.0,7979.0,21.0,10,10);

  // init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001,3.0,1.2,3.8,2.0,16.0,5.0,0.8);
  init_priors(0.002,sigma_zphot_shear[sce],0.001,0.001,sigma_zphot_clustering[sce],0.001,0.001);
  init_survey(survey_designation[sce],nsource_table[sce],nlens_table[sce],area_table[sce]);
  sprintf(arg1,"zdistris/%s",source_zfile[sce]);
  sprintf(arg2,"zdistris/%s",lens_zfile[sce]);
  init_galaxies(arg1,arg2,"gaussian","gaussian",tomo_binning_source[sce],tomo_binning_lens[sce]);
  // init_IA("NLA_HF","GAMA"); 
  init_IA("NLA_z","none"); 
  init_probes("10x2pt");
  init_cmb(cmb_yr[sce]);

  init_feedback(1);

#ifdef ONESAMPLE
  sprintf(arg3,"%s_1sample",survey_designation[sce]);
#else
  sprintf(arg3,"%s",survey_designation[sce]);
#endif

  input_cosmo_params_local p_cosmo;
  input_nuisance_params_local p_sys;

  p_cosmo = (input_cosmo_params_local) {0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.};
  // p_cosmo = (input_cosmo_params_local) {0.3,0.8281663873060578,0.97,-1.,0.,0.05,0.7,0.,0.};
  p_sys = (input_nuisance_params_local) {.source_z_s=sigma_zphot_shear[sce], .lens_z_s=sigma_zphot_clustering[sce], \
                                         .A_ia=0.5, .eta_ia=0.};
  double p_gas[15] = {1.17,0.6,14.,1.,0.03,12.5,1.2,\
                      6.5,0.752,0.,0.,\
                      0.6,14.,0.,0.}; // last 4 for are param_v2 for TEST_CALIB, not used in regular run
  // double p_gas[15] = {1.17702,0.6,13.59369,0.84710,0.0330,12.4479,1.2,\
  //                     6.65445,0.752,-0.10650,0.,\
                         0.6,13.59369,-0.10650,0.}; // HMCODE test T_AGN=1e7.8 K


  for(i=0;i<10;i++){
    p_sys.bias[i] = gbias.b[i];
    p_sys.source_z_bias[i] = 0.;
    p_sys.lens_z_bias[i] = 0.;
    p_sys.shear_m[i] = 0.;
  }
#ifdef TEST_CALIB
  for(i=0;i<15;i++){
#else
  for(i=0;i<11;i++){
#endif
    p_sys.gas[i] = p_gas[i];
  }

  compute_data_vector(arg3,p_cosmo,p_sys);

  // compute_data_vector(arg3,0.3,0.8281663873060578,0.97,-1.,0.,0.05,0.7,0.,0.,\
  //   gbias.b[0],gbias.b[1],gbias.b[2],gbias.b[3],gbias.b[4],\
  //   gbias.b[5],gbias.b[6],gbias.b[7],gbias.b[8],gbias.b[9],\
  //   0.0,0.0,0.0,0.0,0.0,\
  //   0.0,0.0,0.0,0.0,0.0,\
  //   sigma_zphot_shear[sce],\
  //   0.0,0.0,0.0,0.0,0.0,\
  //   0.0,0.0,0.0,0.0,0.0,\
  //   sigma_zphot_clustering[sce],\
  //   0.0,0.0,0.0,0.0,0.0,\
  //   0.0,0.0,0.0,0.0,0.0,\
  //   5.92,-0.47,\
  //   1.17702,0.6,13.59369,0.84710,0.0330,12.4479,1.2,\
  //   6.65445,0.752,-0.10650,0.); // HMCODE test T_AGN=1e7.8 K

  // FILE *f = fopen("pdelta_hm_trans.txt","w");
  // for(double aa=0.6;aa<1.;aa+=0.001){
  //   for(double logk=log(limits.k_min_cH0);logk<log(limits.k_max_cH0);logk+=0.1){
  //     fprintf(f, "%le %le %le\n", aa, exp(logk), Pdelta(exp(logk),aa));
  //   }
  // }
  // fclose(f);

  // FILE *f = fopen("Cl_ss_00_hm_trans.txt","w");
  // // FILE *f = fopen("Cl_ss_00_hf.txt","w");
  // for(double ll=30.;ll<3000.;ll+=100){
  //   fprintf(f, "%le %le\n", ll, C_shear_tomo_sys(ll,0,0));
  // }
  // fclose(f);
 
  // FILE *f = fopen("components_ss_00_hm_trans.txt","w");
  // // FILE *f = fopen("components_ss_00_hf.txt","w");
  // for(double aa=0.1;aa<1.;aa+=0.02){
  //   for(double logk=log(limits.k_max_cH0/100.);logk<log(limits.k_max_cH0);logk+=0.5){
  //     fprintf(f, "%le %le %le\n", aa, exp(logk), Pdelta(exp(logk),aa));
  //   }
  //   // fprintf(f, "%le %le %le %le %le\n", aa, W_kappa(aa,chi(aa),0), A_IA_Joachimi(aa), chi(aa), growfac(aa));
  // }
  // fclose(f);  

  // double kk;
  // FILE *f = fopen("I12_SSC.txt","w");
  // for(double aa=0.2;aa<1.;aa+=0.02){
  //   for(double logk=log(limits.k_min_cH0);logk<log(limits.k_max_cH0);logk+=0.5){
  //     kk = exp(logk);
  //     fprintf(f, "%le %le %le %le %le\n", aa, kk, I12_SSC(kk,aa), I12_SSC_yy(kk,aa),I12_SSC_my(kk,aa));
  //   }
  // }
  // fclose(f);

  return 0;
}


