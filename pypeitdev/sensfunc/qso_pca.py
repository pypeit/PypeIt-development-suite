


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from sklearn import mixture
import pickle
import IPython



def init_pca(filename,wave_grid,redshift, npca):
    # Read in the pickle file from coarse_pca.create_coarse_pca
    # The relevant pieces are the wavelengths (wave_pca_c), the PCA components (pca_comp_c),
    # and the Gaussian mixture model prior (mix_fit)

    loglam = np.log10(wave_grid)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    wave_pca_c, cont_all_c, pca_comp_c, coeffs_c, mean_pca, covar_pca, diff_pca, mix_fit, chi2, dof = pickle.load(open(filename,'rb'))
    num_comp = pca_comp_c.shape[0] # number of PCA components
    # Interpolate PCA components onto wave_grid
    pca_interp = scipy.interpolate.interp1d(wave_pca_c*(1.0 + redshift),pca_comp_c, bounds_error=False, fill_value=0.0, axis=1)
    pca_comp_new = pca_interp(wave_grid)
    # Generate a mixture model for the coefficients prior, what should ngauss be?
    prior = mixture.GaussianMixture(n_components = npca-1).fit(coeffs_c[:, 1:npca])
    # Construct the PCA dict
    pca_dict = {'npca': npca, 'components': pca_comp_new, 'prior': prior, 'coeffs': coeffs_c,
                'z_fid': redshift, 'dloglam': dloglam}
    return pca_dict

def pca_eval(theta,pca_dict):
    # theta_pca[0] is redshift
    # theta_pca[1] is norm
    # theta_pca[2:npca+1] are the PCA dimensionality

    C = pca_dict['components']
    z_fid = pca_dict['z_fid']
    dloglam = pca_dict['dloglam']
    npca = pca_dict['npca']  # Size of the PCA currently being used, original PCA in the dict could be larger
    z_qso = theta[0]
    norm = theta[1]
    A = theta[2:]
    dshift = int(np.round(np.log10((1.0 + z_qso)/(1.0 + z_fid))/dloglam))
    C_now = np.roll(C[:npca,:], dshift, axis=1)
    return norm*np.exp(np.dot(np.append(1.0,A),C_now))

def pca_lnprior(theta,pca_dict):
    gaussian_mixture_model = pca_dict['prior']
    A = theta[2:]
    return gaussian_mixture_model.score_samples(A.reshape(1,-1))