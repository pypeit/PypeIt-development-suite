import numpy as np
from sklearn import mixture

def init_pca(filename,wave_grid,redshift):
    # Read in the pickle file
    # The relevant pieces are the wavelengths (wave_pca_c), the PCA components (pca_comp_c),
    # and the Gaussian mixture model prior (mix_fit)
    wave_pca_c, cont_all_c, pca_comp_c, coeffs_c,
    mean_pca, covar_pca, diff_pca, mix_fit, chi2, dof = pickle.load(open(filename,'rb'))
    num_comp = pca_comp_c.shape[0] # number of PCA components
    # Interpolate PCA components onto wave_grid
    pca_interp = scipy.interpolate.interp1d(wave_pca_c*(1+redshift),pca_comp_c,
                                            bounds_error=False, axis=1)
    pca_comp_new = pca_interp(wave_grid)
    # Construct the PCA dict
    pca_dict = {'mean_spectrum': mean_pca, 'n_components': num_comp,
                'components': pca_comp_new, 'prior': mix_fit}
    return pca_dict

def eval_pca(theta,pca_dict):
    M = pca_dict['mean_spectrum']
    C = pca_dict['components']
    norm = theta[0]
    A = theta[1:]
    return norm*M*np.exp(np.dot(A,C))

def eval_pca_prior(theta,pca_dict):
    gmm = pca_dict['prior']
    A = theta[1:]
    return gmm.score_samples(A.reshape(1,-1))

