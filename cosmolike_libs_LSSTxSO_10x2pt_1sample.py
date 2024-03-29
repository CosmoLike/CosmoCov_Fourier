import sys
import emcee
import ctypes
import os
import numpy as np
import mpi4py

# from mpp_blinding import blind_parameters
# from mpp_blinding import seed as blinding_seed

dirname = os.path.split(__file__)[0]
lib_name = os.path.join(dirname, "like_fourier_10x2pt_fft_1sample.so")
lib=ctypes.cdll.LoadLibrary(lib_name)
double = ctypes.c_double

Double10 = double*10

initcosmo=lib.init_cosmo_runmode
initcosmo.argtypes=[ctypes.c_char_p]

initbins=lib.init_binning_fourier
initbins.argtypes=[ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]

initsurvey=lib.init_survey
initsurvey.argtypes=[ctypes.c_char_p,ctypes.c_double, ctypes.c_double, ctypes.c_double]


initgalaxies=lib.init_galaxies
initgalaxies.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p]

initia=lib.init_IA
initia.argtypes=[ctypes.c_char_p,ctypes.c_char_p]

initfb=lib.init_feedback
initfb.argtypes=[ctypes.c_int]

initpriors=lib.init_priors
initpriors.argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_double]

initprobes=lib.init_probes
initprobes.argtypes=[ctypes.c_char_p]

initcmb=lib.init_cmb
initcmb.argtypes=[ctypes.c_char_p]

initdatainv=lib.init_data_inv_mask
initdatainv.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p]

get_N_tomo_shear = lib.get_N_tomo_shear
get_N_tomo_shear.argtypes = []
get_N_tomo_shear.restype = ctypes.c_int

get_N_tomo_clustering = lib.get_N_tomo_clustering
get_N_tomo_clustering.argtypes = []
get_N_tomo_clustering.restype = ctypes.c_int

get_N_ggl = lib.get_N_ggl
get_N_ggl.argtypes = []
get_N_ggl.restype = ctypes.c_int

get_N_ell = lib.get_N_ell
get_N_ell.argtypes = []
get_N_ell.restype = ctypes.c_int

skip_shearcalib_phz_sampling = lib.skip_shearcalib_phz_sampling
skip_shearcalib_phz_sampling.argtypes = []

# lib.initialize_all_wrapper.restype = ctypes.c_int
# lib.initialize_all_wrapper.argtypes = [
#     ctypes.c_char_p,       # const char * base_dir, 
#     ctypes.c_bool,         # bool auborg_prior, 
#     ctypes.c_bool,         # bool photo_bao_prior, 
#     ctypes.c_bool,         # bool ia_datavector, 
#     ctypes.c_char_p,       # char * ia_model, // must be either "none", "NLA_HF"
#     ctypes.c_char_p,       # char * ia_luminosity_function // should be either "GAMA", "DEEP2"
#     ctypes.c_bool,         # bool modify_shear_priors,
#     Double10,              # double shear_m_mean[10],
#     Double10,              # double shear_m_var[10],
#     ctypes.c_bool,         # bool modify_photoz_priors,
#     Double10,              # double photoz_source_bias_mean[10],
#     double,                # double photoz_source_sigma_mean,
#     Double10,              # double photoz_source_bias_var[10],
#     double,                # double photoz_source_sigma_var,
#     ctypes.c_char_p,       # source_filename
#     ctypes.c_char_p,       # lens_filename
# ]


# initialize_all_wrapper=lib.initialize_all_wrapper

class IterableStruct(ctypes.Structure):
    def names(self):
        out = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                out.append(name)
            else:
                for i in range(length):
                    out.append(name + "_" + str(i))
        return out


    def iter_parameters(self):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                yield name, obj, ptype._length_
            else:
                yield name, obj, 0

    def iter_parameters_filter(self, used):
        for (name, obj, length) in self.iter_parameters():
            if name in used:
                yield name, obj, 0


    def convert_to_vector(self):
        p = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                p.append(obj)
            else:
                for i in range(length):
                    p.append(obj[i])
        return p

    def convert_to_vector_filter(self, used):
        p = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                if name in used:
                    p.append(obj)
            else:
                for i in range(length):
                    if name+'_'+str(i) in used:
                        p.append(obj[i])
        return p



    def read_from_cosmosis(self, block):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                for i in range(ptype._length_):
                    obj[i] = block[self.section_name, name+"_"+str(i)]
            else:
                setattr(self, name, block[self.section_name, name])



    def print_struct(self):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                for i in range(ptype._length_):
                    print("%s[%d] = %f"%(name, i, obj[i]))
            else:
                print("%s = %f"%(name, obj))


    def number_of_doubles(self):
        n=0
        for name, ptype in self._fields_:
            if hasattr(ptype, "_length_"):
                n += ptype._length_
            else:
                n += 1
        return n

    def set_from_vector(self, p):
        i=0
        j=0
        while i<len(p):
            name,ptype = self._fields_[j]
            j+=1
            if ptype == double:
                setattr(self, name, p[i])
                i+=1
            else:
                x = getattr(self, name)
                assert x._type_==double
                for k in range(x._length_):
                    x[k] = p[i]
                    i+=1


class InputCosmologyParams(IterableStruct):
    section_name = "cosmological_parameters"
    _fields_ = [
        ("omega_m", double),
        ("sigma_8", double),
        ("n_s", double),
        ("w0", double),
        ("wa", double),
        ("omega_b", double),
        ("h0", double),
        ("MGSigma", double),
        ("MGmu", double),
    ]

    @classmethod
    def fiducial(cls):
        c = cls()
        c.omega_m = 0.3156
        c.sigma_8 = 0.831
        c.n_s = 0.9645
        c.w0 = -1.0
        c.wa = 0.0
        c.omega_b = 0.0491685
        c.h0 = 0.6727
        c.MGSigma = 0.0
        c.MGmu = 0.0
        return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.omega_m = 0.02
        c.sigma_8 = 0.02
        c.n_s = 0.01
        c.w0 = 0.1
        c.wa = 0.1
        c.omega_b = 0.001
        c.h0 = 0.01
        c.MGSigma = 0.1
        c.MGmu = 0.1     
        return c



class InputNuisanceParams(IterableStruct):
    section_name = "nuisance_parameters"
    _fields_ = [
        ("bias", double*10),
        ("source_z_bias", double*10),
        ("source_z_s", double),
        ("lens_z_bias", double*10),
        ("lens_z_s", double),
        ("shear_m", double*10),
        ("A_ia", double),
        ("eta_ia", double),
        ("gas", double*11)
    ]
    @classmethod
    def fiducial(cls):
        c = cls()
        c.bias[:] = [1.125013,1.308513,1.433975,1.558521,1.693127,1.844866,2.026988,2.265804,2.633997,3.6] # Y6; fid value can't be too close to prior boundary, sample will go over!
        c.source_z_bias[:] = np.repeat(0.0, 10)
        c.source_z_s = 0.05
        c.lens_z_bias[:] = np.repeat(0.0, 10)
        c.lens_z_s = 0.03
        c.shear_m[:] = np.repeat(0.0, 10)
        c.A_ia = 0.5
        c.eta_ia = 0.
        c.gas[:] = [1.17,0.6,14.,1.,0.03,12.5,1.2,\
                    6.5,0.752,0.,0.]
        return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.bias[:] = np.repeat(0.2, 10)
        c.source_z_bias[:] = np.repeat(0.01, 10)
        c.source_z_s = 0.005
        c.lens_z_bias[:] = np.repeat(0.01, 10)
        c.lens_z_s = 0.005
        c.shear_m[:] = np.repeat(0.01, 10)
        c.A_ia = 0.1
        c.eta_ia = 0.1
        c.gas[:] = [0.01,0.01,0.05,0.05,0.001,0.5,0.05,\
                    0.05,0.01,0.05,0.05]
        return c

class LikelihoodFunctionWrapper(object):
    def __init__(self, varied_parameters):
        self.varied_parameters = varied_parameters


    def fill_varied(self, icp, inp, x):
        assert len(x) == len(self.varied_parameters), "Wrong number of parameters"
        i = 0
        for s in [icp, inp]:
            for name, obj, length in s.iter_parameters():
                if length==0:
                    if name in self.varied_parameters:
                        setattr(s, name, x[i])        
                        i+=1
                else:
                    for j in range(length):
                        name_i = name + "_" + str(j)
                        if name_i in self.varied_parameters:
                            obj[j] = x[i]
                            i+=1

    def __call__(self, x):
        icp = InputCosmologyParams.fiducial()
        inp = InputNuisanceParams.fiducial()
        self.fill_varied(icp, inp, x)
        #icp.print_struct()
        #inp.print_struct()
        #print
        like = lib.log_multi_like(icp, inp)
        return like

lib.log_multi_like.argtypes = [InputCosmologyParams, InputNuisanceParams]
lib.log_multi_like.restype = double
log_multi_like = lib.log_multi_like

lib.compute_data_vector.argtypes = [ctypes.c_char_p, InputCosmologyParams, InputNuisanceParams]
compute_data_vector = lib.compute_data_vector

def sample_cosmology_only(MG = False, w0wa = False):

    varied_parameters = ['omega_m']
    varied_parameters.append('sigma_8')
    varied_parameters.append('n_s')

    if(w0wa):
        varied_parameters.append('w0')
        varied_parameters.append('wa')

    varied_parameters.append('omega_b')
    varied_parameters.append('h0')

    if(MG):
        varied_parameters.append('MGSigma')
        varied_parameters.append('MGmu')

    return varied_parameters

def sample_cosmology_10x2_fixhalo(tomo_N_shear,tomo_N_lens,MG = False, w0wa=False, cov_modified=False):
    varied_parameters = sample_cosmology_only(MG,w0wa)
    varied_parameters += ['bias_%d'%i for i in range(tomo_N_lens)]
    if cov_modified is False:
        varied_parameters += ['source_z_bias_%d'%i for i in range(tomo_N_shear)]
        varied_parameters.append('source_z_s')
        # varied_parameters += ['lens_z_bias_%d'%i for i in range(tomo_N_lens)]
        # varied_parameters.append('lens_z_s')
        varied_parameters += ['shear_m_%d'%i for i in range(tomo_N_shear)]

    varied_parameters.append('A_ia')
    varied_parameters.append('eta_ia')

    return varied_parameters

def sample_cosmology_10x2_fix_eps12Gamma(tomo_N_shear,tomo_N_lens,MG = False, w0wa=False, cov_modified=False):
    varied_parameters = sample_cosmology_only(MG,w0wa)
    varied_parameters += ['bias_%d'%i for i in range(tomo_N_lens)]
    if cov_modified is False:
        varied_parameters += ['source_z_bias_%d'%i for i in range(tomo_N_shear)]
        varied_parameters.append('source_z_s')
        # varied_parameters += ['lens_z_bias_%d'%i for i in range(tomo_N_lens)]
        # varied_parameters.append('lens_z_s')
        varied_parameters += ['shear_m_%d'%i for i in range(tomo_N_shear)]

    varied_parameters.append('A_ia')
    varied_parameters.append('eta_ia')

    i_gas = [1,2,3, 7,8] # select gas parameters to vary
    varied_parameters += ['gas_%d'%i for i in i_gas]
    return varied_parameters

def sample_cosmology_10x2_allsys(tomo_N_shear,tomo_N_lens,MG = False, w0wa=False, cov_modified=False):
    varied_parameters = sample_cosmology_only(MG,w0wa)
    varied_parameters += ['bias_%d'%i for i in range(tomo_N_lens)]
    if cov_modified is False:
        varied_parameters += ['source_z_bias_%d'%i for i in range(tomo_N_shear)]
        varied_parameters.append('source_z_s')
        # varied_parameters += ['lens_z_bias_%d'%i for i in range(tomo_N_lens)]
        # varied_parameters.append('lens_z_s')
        varied_parameters += ['shear_m_%d'%i for i in range(tomo_N_shear)]

    varied_parameters.append('A_ia')
    varied_parameters.append('eta_ia')

    i_gas = [0,1,2,3, 7,8,9,10] # select gas parameters to vary
    varied_parameters += ['gas_%d'%i for i in i_gas]
    return varied_parameters

def sample_cosmology_3x2_allsys(tomo_N_shear,tomo_N_lens,MG = False, w0wa=False, cov_modified=False):
    varied_parameters = sample_cosmology_only(MG,w0wa)
    varied_parameters += ['bias_%d'%i for i in range(tomo_N_lens)]
    if cov_modified is False:
        varied_parameters += ['source_z_bias_%d'%i for i in range(tomo_N_shear)]
        varied_parameters.append('source_z_s')
        # varied_parameters += ['lens_z_bias_%d'%i for i in range(tomo_N_lens)]
        # varied_parameters.append('lens_z_s')
        varied_parameters += ['shear_m_%d'%i for i in range(tomo_N_shear)]

    varied_parameters.append('A_ia')
    varied_parameters.append('eta_ia')

    i_gas = [1,2,9,10] # select gas parameters to vary
    varied_parameters += ['gas_%d'%i for i in i_gas]
    return varied_parameters

# def sample_cosmology_3x2_hfsys(tomo_N_shear,tomo_N_lens,MG = False, w0wa=False, cov_modified=False):
#     varied_parameters = sample_cosmology_only(MG,w0wa)
#     varied_parameters += ['bias_%d'%i for i in range(tomo_N_lens)]
#     if cov_modified is False:
#         varied_parameters += ['source_z_bias_%d'%i for i in range(tomo_N_shear)]
#         varied_parameters.append('source_z_s')
#         varied_parameters += ['lens_z_bias_%d'%i for i in range(tomo_N_lens)]
#         varied_parameters.append('lens_z_s')
#         varied_parameters += ['shear_m_%d'%i for i in range(tomo_N_shear)]

#     varied_parameters.append('A_ia')
#     varied_parameters.append('eta_ia')

#     # i_gas = [1,2,9,10] # select gas parameters to vary
#     # varied_parameters += ['gas_%d'%i for i in i_gas]
#     return varied_parameters

# def sample_cosmology_ss_hfsys(tomo_N_shear,MG = False, w0wa=False, cov_modified=False):
#     varied_parameters = sample_cosmology_only(MG,w0wa)
#     # varied_parameters += ['bias_%d'%i for i in range(tomo_N_lens)]
#     if cov_modified is False:
#         varied_parameters += ['source_z_bias_%d'%i for i in range(tomo_N_shear)]
#         varied_parameters.append('source_z_s')
#         # varied_parameters += ['lens_z_bias_%d'%i for i in range(tomo_N_lens)]
#         # varied_parameters.append('lens_z_s')
#         varied_parameters += ['shear_m_%d'%i for i in range(tomo_N_shear)]

#     varied_parameters.append('A_ia')
#     varied_parameters.append('eta_ia')

#     # i_gas = [1,2,9,10] # select gas parameters to vary
#     # varied_parameters += ['gas_%d'%i for i in i_gas]
#     return varied_parameters


def sample_main(varied_parameters,sigma_z_shear,sigma_z_clustering, iterations, nwalker, nthreads, filename, blind=False, pool=None):
    print(varied_parameters)

    likelihood = LikelihoodFunctionWrapper(varied_parameters)
    starting_point = InputCosmologyParams.fiducial().convert_to_vector_filter(varied_parameters)
    
    #changing the center of the 'starting sphere' of the MCMC, according to the fiducial input parameter 
    new=InputNuisanceParams().fiducial()
    setattr(new,'source_z_s',sigma_z_shear)
    setattr(new,'lens_z_s',sigma_z_clustering) 
    starting_point += new.convert_to_vector_filter(varied_parameters)
    #starting_point += InputCosmologyParams.fiducial().convert_to_vector_filter(varied_parameters)

    std = InputCosmologyParams.fiducial_sigma().convert_to_vector_filter(varied_parameters)
    std += InputNuisanceParams().fiducial_sigma().convert_to_vector_filter(varied_parameters)

    p0 = emcee.utils.sample_ball(starting_point, std, size=nwalker)

    ndim = len(starting_point)
    print("ndim = %d"%(ndim))
    print("start = ", starting_point)
    print("std = ", std)


    # if pool is not None:
    #     if not pool.is_master():
    #         pool.wait()
    #         sys.exit(0)


    sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood,threads=nthreads,pool=pool)

#    sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood, pool=pool)

    f = open(filename, 'w')

    #write header here
    f.write('# ' + '    '.join(varied_parameters)+" log_like\n")
    f.write('#blind=%s\n'%blind)
    if blind:
        f.write('#blinding_seed=%d\n'%blinding_seed)

    for (p, loglike, state) in sampler.sample(p0,iterations=iterations):
        for row,logl in zip(p,loglike):
            if blind:
                row = blind_parameters(varied_parameters, row)
            p_text = '  '.join("%.5e"%r for r in row)
            f.write('%s %e\n' % (p_text,logl))
        f.flush()
    f.close()
    

    # for (p, loglike, state) in sampler.sample(p0,iterations=iterations):
    #     for row in p:
    #         if blind:
    #             row = blind_parameters(varied_parameters, row)
    #         p_text = '  '.join(str(r) for r in row)
    #         print ('%s %e\n' % (p_text,loglike))
    #         f.write('%s %e\n' % (p_text,loglike))
    #     f.flush()
    # f.close()

