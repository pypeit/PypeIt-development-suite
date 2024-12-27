import numpy as np
import scipy
from sklearn import mixture

def init_pca(filename,wave_grid,redshift):
    # Read in the pickle file from coarse_pca.create_coarse_pca
    # The relevant pieces are the wavelengths (wave_pca_c), the PCA components (pca_comp_c),
    # and the Gaussian mixture model prior (mix_fit)
    wave_pca_c, cont_all_c, pca_comp_c, coeffs_c,
    mean_pca, covar_pca, diff_pca, mix_fit, chi2, dof = pickle.load(open(filename,'rb'))
    num_comp = pca_comp_c.shape[0] # number of PCA components
    # Interpolate PCA components onto wave_grid
    pca_interp = scipy.interpolate.interp1d(wave_pca_c*(1+redshift),pca_comp_c,
                                            bounds_error=False, axis=1)
    nord = wave_grid.shape[1]
    pca_comp_new = np.zeros(wave_grid.shape)
    for ii in range(nord):
        pca_comp_new[ii] = pca_interp(wave_grid[:,ii])
    # Construct the PCA dict
    pca_dict = {'n_components': num_comp, 'components': pca_comp_new,
            'prior': mix_fit, 'coeffs': coeffs_c}
    return pca_dict

def eval_pca(theta,pca_dict):
    C = pca_dict['components']
    norm = theta[0]
    A = theta[1:]
    return norm*np.exp(np.dot(np.append(1.0,A),C))

def eval_pca_prior(theta,pca_dict):
    gmm = pca_dict['prior']
    A = theta[1:]
    return gmm.score_samples(A.reshape(1,-1))

