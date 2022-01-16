cfastcov_dir := cfastcov/
cfastcov := $(cfastcov_dir)twobessel.c $(cfastcov_dir)utils.c $(cfastcov_dir)utils_complex.c
opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-std=gnu99 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass
opt_ocelote := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
opt_puma := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c

home_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_home)

home_cov_102pt:
	gcc compute_covariances_fourier_10x2pt.c -o ./compute_covariances_fourier_10x2pt $(opt_home)


home_cov_1sample:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier_1sample $(opt_home) -DONESAMPLE

home_datav_102pt:
	gcc like_fourier_10x2pt.c -o ./like_fourier_10x2pt $(opt_home)
# 	gcc like_fourier_10x2pt.c -o ./like_fourier_10x2pt_1sample $(opt_home) -DONESAMPLE


ocelote_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_ocelote)

ocelote_cov_1sample:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier_1sample $(opt_ocelote) -DONESAMPLE

###### Puma

puma_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_puma)
