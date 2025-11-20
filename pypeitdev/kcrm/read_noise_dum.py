
import numpy as np

sigma_rn = 3.9
slit_factor = 4.0
binning_factor = 4.0
sky_rate_1x1 = 80.0/slit_factor/binning_factor/300.0 # sky rate in e-/s/pix 1x1 binning
sky_rate_2x2 = 80.0/slit_factor/300.0 # sky rate in e-/s/pix 2x2 binning

t_total = 3600.0

# consider a 5 exposure sequence
t1 = 900.0
nexp1 = t_total/t1
sigma1 = np.sqrt(nexp1)*np.sqrt(sky_rate_1x1*t1 + sigma_rn**2)

# consider a 6 exposure sequence with 600s
t2 = 300.0
nexp2 = t_total/t2
sigma2 = np.sqrt(nexp2)*np.sqrt(sky_rate_2x2*t2 + sigma_rn**2)

# Consider an upscaling of the strategy 1
t_total_alt=1.17*t_total
nalt = t_total_alt/t1
sigma_alt = np.sqrt(nalt)*np.sqrt(sky_rate_1x1*t1 + sigma_rn**2)

obj_rate_2 = 0.3*sky_rate_2x2


f_total = obj_rate_2*t_total # total flux in a pixel of size 2x2
f_total_alt = obj_rate_2*t_total_alt

SNR1 = f_total/(np.sqrt(binning_factor)*sigma1) # 4 times more pixels for the 1x1 binned data
SNR2 = f_total/sigma2
SNR_alt = f_total_alt/(np.sqrt(binning_factor)*sigma_alt)




