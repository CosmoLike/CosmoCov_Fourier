import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import pandas as pd

class Mask:
    def __init__(self, path = None):
        if path:
            self.set_mask(path)
        else:
            print("Mask object must be initialized with file path")
    
    def set_mask(self, path):
        mask_file= np.loadtxt(path)
        self.mask = np.int_(mask_file[:,1])
        print("mask set.")

class RawCov:
    def __init__(self, path = None, isGaussian=False):
        if path:
            self.covmat = self.convert_cov(path, isGaussian)
        else:
            print("RawCov object must be initialized with file path")
            
    def convert_cov(self, filename, isGaussian):
        print("reading cov: %s"%(filename))
        data = pd.read_csv(filename, sep=' ', header=None, usecols=[0,1,8,9])
        ndata = int(np.max(data.iloc[:,0]))+1
        print("row counts: %d"%data.shape[0])
        print("size: %d x %d"%(ndata, ndata))
        cov_g = np.zeros((ndata,ndata))
        cov_ng = np.zeros((ndata,ndata))
    #     for i in range(0,data.shape[0]):
    #         x = int(data.iat[i,0])
    #         y = int(data.iat[i,1])
    #         cov_g[x,y] = cov_g[y,x] = data.iat[i,2]
    #         cov_ng[x,y] = cov_ng[y,x] = cov_g[x,y] + data.iat[i,3]
        x = np.int_(np.array(data.iloc[:,0]))
        y = np.int_(np.array(data.iloc[:,1]))
        if isGaussian:
            cov_g[x, y] = np.array(data.iloc[:,2])
            cov_g[y, x] = cov_g[x, y]
        cov_ng[x, y] = np.array(data.iloc[:,2]) + np.array(data.iloc[:,3])
        cov_ng[y, x] = cov_ng[x, y]
        np.savetxt(filename+'_comp', np.c_[np.repeat(np.arange(ndata),ndata), np.tile(np.arange(ndata),ndata), cov_ng.flatten()], fmt="%d %d %le")
        
        if isGaussian:
            np.savetxt(filename+'_gaussian_comp', np.c_[np.repeat(np.arange(ndata),ndata), np.tile(np.arange(ndata),ndata), cov_g.flatten()], fmt="%d %d %le")
        print('Compressed covariance file saved at %s'%(filename+'_comp'))

        return cov_g if isGaussian else cov_ng

class Cov:
    
    def __init__(self, path = None, matrix = None):
        if path:
            self.read_mat(path)
        else:
            self.covmat = matrix
            self.ndata = self.covmat.shape[0]
            print("covmat assigned, dimension: %d x %d"%(self.ndata, self.ndata))
            
    def read_mat(self, filename):
        print("reading mat: %s"%(filename))
        data = pd.read_csv(filename, sep=' ', header=None, usecols=[0,1,2])
        self.ndata = int(np.max(data.iloc[:,0]))+1
        print("row counts: %d"%data.shape[0])
        print("size: %d x %d"%(self.ndata, self.ndata))
        self.covmat = np.array(data.iloc[:,2]).reshape((self.ndata, self.ndata))
        
    def get_corr_mat(self):
        corr = np.zeros(self.covmat.shape)
        for i in range(self.ndata):
            for j in range(self.ndata):
                corr[i][j] = self.covmat[i][j] / np.sqrt(self.covmat[i][i]*self.covmat[j][j])
        return corr
    
    def plot_mat(self, figsize=(15,15), cmap="seismic", vmin=-1, vmax=1, savefig=None, dpi=300):
        plt.figure(figsize=figsize)
        im = plt.imshow(self.covmat, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, orientation='vertical')
        if savefig is not None:
            plt.savefig(savefig,dpi=dpi)
        plt.show()

    def get_invcov_masked(self, mask=None, inds=None):
        if mask is None:
            cov_partmasked = self.covmat
        else:
            maskinds = (mask == 1)
            cov_diag= np.diag(np.diag(self.covmat))
            cov_offdiag = self.covmat - cov_diag
            cov_offdiag[maskinds,:][:,maskinds] = 0.
            cov_partmasked = cov_offdiag + cov_diag
        
        if inds is None:
            invc = LA.inv(cov_partmasked)
        else:
            invc = LA.inv(cov_partmasked[inds,:][:,inds])
            
        return invc
    
    def get_det_eig_cov_masked(self, mask=None, inds=None):
        fullmask = np.ones(self.ndata, dtype=np.int)
        masked = False
        if mask is not None:
            fullmask *= mask
            masked = True
        
        if inds is not None:
            tmpmask = np.zeros(self.ndata, dtype=np.int)
            tmpmask[inds] = 1
            fullmask *= tmpmask
            masked = True
        
        fullinds = (fullmask==1)
        detc = LA.det(self.covmat[fullinds,:][:,fullinds])
        # if not masked:
        #     eigens = LA.eigvals(self.covmat)
        # else:
        eigens = LA.eigvals(self.covmat[fullinds,:][:,fullinds])
            
        return detc, eigens

    def test_eigens(self, mask=None, inds=None):
        _, eigens = self.get_det_eig_cov_masked(mask, inds)
        a = np.sort(eigens)
        print('sorted eigenvalues: ', a[:])
        if a[0] <= 0:
            raise ValueError('Non-positive eigenvalues!')
        print('Eigenvalues look good!')
    
    def get_snr(self, dvfile, mask=None, inds=None):
        dv = np.loadtxt(dvfile)[:,1]
        invc = self.get_invcov_masked(mask, inds)
        return np.sqrt(invc.dot(dv).dot(dv))

def cov_3col_out(mat, outfile):
    n = mat.shape[0]
    narr = np.arange(n)
    col1 = np.repeat(narr, n)
    col2 = np.tile(narr, n)
    np.savetxt(outfile, np.c_[col1, col2, mat.flatten()], fmt='%d %d %le')
    print('cov saved to %s!'%(outfile))