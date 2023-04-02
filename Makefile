cfastcov_dir := cfastcov/
cfastcov := $(cfastcov_dir)twobessel.c $(cfastcov_dir)utils.c $(cfastcov_dir)utils_complex.c
opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -O3 \
-std=gnu99 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass
opt_ocelote := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lfftw3 -lgsl -lgslcblas -lm -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
opt_puma := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib \
-lfftw3 -lgsl -lgslcblas -lm -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c

home_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_home)

home_cov_102pt:
	gcc compute_covariances_fourier_10x2pt.c -o ./compute_covariances_fourier_10x2pt $(opt_home)
home_cov_102pt_slow:
	gcc compute_covariances_fourier_10x2pt.c -o ./compute_covariances_fourier_10x2pt_slow $(opt_home) -DSLOW

home_cov_102pt_fft:
	gcc compute_covariances_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -o ./compute_covariances_fourier_10x2pt_fft $(opt_home)

home_cov_102pt_fft_1sample:
	gcc compute_covariances_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DONESAMPLE -o ./compute_covariances_fourier_10x2pt_fft_1sample $(opt_home)

home_cov_102pt_fft_test:
	gcc compute_covariances_fourier_10x2pt_test.c $(cfftlog) -DRUN_FFT -o ./compute_covariances_fourier_10x2pt_fft_test $(opt_home)

home_cov_3x2pt_fft_hf:
	gcc compute_covariances_fourier_3x2pt_halofit.c $(cfftlog) -DRUN_FFT -o ./compute_covariances_fourier_3x2pt_fft_hf $(opt_home)


home_cov_1sample:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier_1sample $(opt_home) -DONESAMPLE

home_datav_102pt:
	gcc like_fourier_10x2pt.c $(cfftlog) -o ./like_fourier_10x2pt $(opt_home)
# 	gcc like_fourier_10x2pt.c -o ./like_fourier_10x2pt_1sample $(opt_home) -DONESAMPLE

home_datav_102pt_fft:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -o ./like_fourier_10x2pt_fft $(opt_home)

# use Ludlow16 mass-concentration for datav
home_datav_102pt_fft_Ludlow16:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DLUDLOW16 -o ./like_fourier_10x2pt_fft_Ludlow $(opt_home)
home_datav_102pt_fft_Ludlow16fit:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DLUDLOW16_FIT -o ./like_fourier_10x2pt_fft_Ludlow_fit $(opt_home)

# use Castro HMF
home_datav_102pt_fft_Castro:
	make home_datav_102pt_fft_Castro_ROCKSTAR
	make home_datav_102pt_fft_Castro_AHF
	make home_datav_102pt_fft_Castro_SUBFIND
	make home_datav_102pt_fft_Castro_VELOCI

home_datav_102pt_fft_Castro_ROCKSTAR:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DHMF_CASTRO_ROCKSTAR -o ./like_fourier_10x2pt_fft_HMF_CASTRO_ROCKSTAR $(opt_home)
home_datav_102pt_fft_Castro_AHF:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DHMF_CASTRO_AHF -o ./like_fourier_10x2pt_fft_HMF_CASTRO_AHF $(opt_home)
home_datav_102pt_fft_Castro_SUBFIND:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DHMF_CASTRO_SUBFIND -o ./like_fourier_10x2pt_fft_HMF_CASTRO_SUBFIND $(opt_home)
home_datav_102pt_fft_Castro_VELOCI:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DHMF_CASTRO_VELOCI -o ./like_fourier_10x2pt_fft_HMF_CASTRO_VELOCI $(opt_home)



home_datav_102pt_fft_1sample:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DONESAMPLE -o ./like_fourier_10x2pt_fft_1sample $(opt_home)

home_datav_102pt_fft_notab:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DNOTAB -o ./like_fourier_10x2pt_fft_notab $(opt_home)

home_datav_102pt_hf:
	gcc like_fourier_10x2pt_halofit.c $(cfftlog) -DRUN_FFT -o ./like_fourier_10x2pt_hf $(opt_home)

home_102pt_fft_shared:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT  -shared -o like_fourier_10x2pt_fft.so -fPIC $(opt_home)

home_102pt_fft_shared_1sample:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DONESAMPLE  -shared -o like_fourier_10x2pt_fft_1sample.so -fPIC $(opt_home)

home_102pt_fft_shared_calib:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DTEST_CALIB  -shared -o like_fourier_10x2pt_fft_nocalib.so -fPIC $(opt_home)

home_102pt_fft_shared_halfcalib:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DTEST_CALIB -DHALFCALIB  -shared -o like_fourier_10x2pt_fft_halfcalib.so -fPIC $(opt_home)

home_tests:
	gcc like_fourier_10x2pt_tests.c $(cfftlog) -o ./like_tests $(opt_home)

ocelote_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_ocelote)

ocelote_cov_1sample:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier_1sample $(opt_ocelote) -DONESAMPLE

###### Puma

puma_102pt_fft_shared:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT  -shared -o like_fourier_10x2pt_fft.so -fPIC $(opt_puma)

puma_102pt_fft_shared_1sample:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DONESAMPLE  -shared -o like_fourier_10x2pt_fft_1sample.so -fPIC $(opt_puma)

puma_102pt_fft_shared_calib:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DTEST_CALIB  -shared -o like_fourier_10x2pt_fft_nocalib.so -fPIC $(opt_puma)

puma_102pt_fft_shared_halfcalib:
	gcc like_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DTEST_CALIB -DHALFCALIB  -shared -o like_fourier_10x2pt_fft_halfcalib.so -fPIC $(opt_puma)

puma_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_puma)

puma_cov_102pt:
	gcc compute_covariances_fourier_10x2pt.c -o ./compute_covariances_fourier_10x2pt $(opt_puma)

puma_cov_102pt_fft:
	gcc compute_covariances_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -o ./compute_covariances_fourier_10x2pt_fft $(opt_puma)

puma_cov_102pt_fft_1sample:
	gcc compute_covariances_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DONESAMPLE -o ./compute_covariances_fourier_10x2pt_fft_1sample $(opt_puma)

puma_cov_102pt_ssc:
	gcc compute_covariances_fourier_10x2pt.c $(cfftlog) -DRUN_FFT -DSSCONLY -o ./compute_covariances_fourier_10x2pt_fft_ssc $(opt_puma)


puma_cov_102pt_slow:
	gcc compute_covariances_fourier_10x2pt.c -o ./compute_covariances_fourier_10x2pt_slow $(opt_puma) -DSLOW

puma_cov_3x2pt_fft_hf:
	gcc compute_covariances_fourier_3x2pt_halofit.c $(cfftlog) -DRUN_FFT -o ./compute_covariances_fourier_3x2pt_fft_hf $(opt_puma)
