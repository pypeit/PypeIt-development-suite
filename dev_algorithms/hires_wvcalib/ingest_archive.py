

import os
import numpy as np
from astropy.table import Table
from pkg_resources import resource_filename
from matplotlib import pyplot as plt
from pypeit.core.wavecal import templates
from pypeit.core.wavecal import wvutils

# Read template file
templ_table_file = os.path.join(
    resource_filename('pypeit', 'data'), 'arc_lines',
    'hires', 'hires_templ.dat')
tbl = Table.read(templ_table_file, format='ascii')
nrows = len(tbl)

order_min = tbl['IOrder'].min()
order_max = 118

order_vec = np.arange(order_min, order_max +1, 1)
norders = order_vec.size

# Subset of orders in every file. Populated indicates whether a given order is populated
lambda_cen = np.zeros((norders, nrows))
ech_angle = np.zeros((norders, nrows))
populated = np.zeros((norders, nrows), dtype=bool)
XDISP_is_red = np.zeros((norders, nrows), dtype=bool)
binspec = np.zeros((norders, nrows), dtype=int)
det = np.zeros((norders, nrows), dtype=int)
xd_angle = np.zeros((norders, nrows))

bluest_order = np.zeros(nrows, dtype=int)
xd_angle_file = np.zeros(nrows)
ech_angle_file = np.zeros(nrows)
det_file = np.zeros(nrows)
XDISP_is_red_file = np.zeros(nrows, dtype=bool)

for irow in np.arange(nrows):
    this_order_vec, this_wave, this_arc = templates.xidl_hires(
        os.path.join(os.getenv('HIRES_CALIBS'), 'ARCS', tbl[irow]['Name']), specbin=tbl[irow]['Rbin'])
    if irow == 0:
        nspec = this_wave.shape[1]
        wave = np.zeros((norders, nrows, nspec))
        arcspec = np.zeros((norders, nrows, nspec))
    else:
        assert this_wave.shape[1] == nspec
    indx = this_order_vec - order_min
    populated[indx, irow] = True
    ech_angle[indx, irow] = tbl[irow]['ECH']
    xd_angle[indx, irow] = tbl[irow]['XDAng']
    XDISP_is_red[indx, irow] = tbl[irow]['XDISP'] == 'RED'
    binspec[indx, irow] =  tbl[irow]['Rbin']
    det[indx, irow] =  tbl[irow]['Chip']

    wave[indx, irow, :] = this_wave
    arcspec[indx, irow, :] = this_arc
    lambda_cen[indx, irow] = np.median(this_wave, axis=1)
    # file specific
    bluest_order[irow] = this_order_vec[-1]
    ech_angle_file[irow] = tbl[irow]['ECH']
    xd_angle_file[irow] = tbl[irow]['XDAng']
    det_file[irow] = tbl[irow]['Chip']
    XDISP_is_red_file[irow] = tbl[irow]['XDISP'] == 'RED'

#all_dlam = []
#all_lam = []
#all_orders = []
for indx, iorder in enumerate(order_vec):
    if np.any(populated[indx, :]):
        this_ech = ech_angle[indx, populated[indx, :]]
        this_xd_angle = xd_angle[indx, populated[indx, :]]
        this_lambda_cen = lambda_cen[indx, populated[indx, :]]
        this_wave = wave[indx, populated[indx, :], :]
        dlam = []
        lam = []
        xlam = []
        for iwave in this_wave:
            dlam += list(wvutils.get_delta_wave(iwave, np.ones_like(iwave,dtype=bool)))
            lam += iwave.tolist()
            #xlam += ((np.array(iwave) - np.array(iwave).min())/(np.array(iwave).max() - np.array(iwave).min())).tolist()
        plt.plot(xlam, 3.0e5*np.array(dlam)/np.array(lam), '.') #, label=f'order={iorder}')
        plt.show()
        #plt.legend()

 #       all_dlam += dlam
 #       all_lam += lam
 #       all_orders += [iorder]*len(lam)


#plt.plot(all_lam, 3.0e5*np.array(all_dlam)/np.array(all_lam), '.')
#plt.legend()
#plt.show()

# Plot the central wavelength vs echelle angle order by order
for indx, iorder in enumerate(order_vec):
    if np.any(populated[indx, :]):
        this_ech = ech_angle[indx, populated[indx, :]]
        this_xd_angle = xd_angle[indx, populated[indx, :]]
        this_lambda_cen = lambda_cen[indx, populated[indx, :]]
        plt.plot(this_ech, this_lambda_cen, 'k.', label=f'order={iorder}')
        plt.legend()
        plt.show()


for xdisp in ['UV', 'RED']:
    for idet in [1,2,3]:
        indx = (XDISP_is_red_file == (xdisp == 'RED')) & (det_file == idet)
        plt.plot(xd_angle_file[indx], bluest_order[indx], 'k.', label=f'XDISP={xdisp}, det={idet}')
        plt.legend()
        plt.show()





