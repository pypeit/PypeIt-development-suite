import numpy as np
import os

from astropy.table import Table
from astropy.io import ascii,fits

import deimos_tools
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

import scipy.ndimage as scipynd
from scipy.optimize import curve_fit

import linetools.utils

# NEED TO GENERALIZE THIS
DEIMOS_DROPBOX = '/Users/mgeha/Dropbox/DEIMOS/'


#######################################################        
# SLIT TABLE
#######################################################
def mk_slit_table(nslits,hdu):


    # ADD SLIT SPECIFC
    cols = [deimos_tools.filled_column('rname','                       ',nslits),
            deimos_tools.filled_column('bname','                       ',nslits),
            deimos_tools.filled_column('rspat',-1.,nslits),
            deimos_tools.filled_column('bspat',-1.,nslits),
            deimos_tools.filled_column('rdet',-1,nslits),
            deimos_tools.filled_column('bdet',-1,nslits),
            deimos_tools.filled_column('rSN',-1.,nslits),
            deimos_tools.filled_column('bSN',-1.,nslits),
            deimos_tools.filled_column('xpos',-1.,nslits),
            deimos_tools.filled_column('ypos',-1.,nslits),
            deimos_tools.filled_column('fit_slope',-1.,nslits),
            deimos_tools.filled_column('fit_b',-1.,nslits),
            deimos_tools.filled_column('fit_los',-1.,nslits),

            deimos_tools.filled_column('objra',-1.,nslits),
            deimos_tools.filled_column('objdec',-1.,nslits),
            deimos_tools.filled_column('slittyp',' ',nslits),
            deimos_tools.filled_column('slitname','                  ',nslits),
            deimos_tools.filled_column('maskdef_id','                ',nslits),

            deimos_tools.filled_column('dslitid',-1,nslits),
            deimos_tools.filled_column('desid',-1,nslits),
            deimos_tools.filled_column('slitx1',-1.,nslits),
            deimos_tools.filled_column('slity1',-1.,nslits),


            deimos_tools.filled_column('rms_arc_r',-1.,nslits),
            deimos_tools.filled_column('rms_arc_b',-1.,nslits),
            deimos_tools.filled_column('rms_sky',-1.,nslits),
           
            deimos_tools.filled_column('telluric_h2o',-1.,nslits),
            deimos_tools.filled_column('telluric_o2',-1.,nslits),
            deimos_tools.filled_column('telluric_h2o_err',-1.,nslits),
            deimos_tools.filled_column('telluric_o2_err',-1.,nslits),
            deimos_tools.filled_column('telluric_w',-1.,nslits),
            deimos_tools.filled_column('telluric_chi2',-1.,nslits),  
            
            
            deimos_tools.filled_column('fmodel_v',-1.,nslits),
            deimos_tools.filled_column('fmodel_w',-1.,nslits),
            deimos_tools.filled_column('fmodel_v_err',-1.,nslits),
            deimos_tools.filled_column('fmodel_w_err',-1.,nslits),
            deimos_tools.filled_column('fmodel_f_acc',-1.,nslits),    # acceptence, should be greater than 0.3
            deimos_tools.filled_column('fmodel_pfile',' ',nslits)

           ]

    slits = Table(cols)
    return slits


#######################################################
# CALC SN
#######################################################
def calc_rb_SN(r1,b1, hdu):

    # TRY SPECTRUM AND CALCULATE SN
    try:
        rSN = np.median(hdu[r1].data['OPT_COUNTS'] * np.sqrt(hdu[r1].data['OPT_COUNTS_IVAR']))
        bSN = np.median(hdu[b1].data['OPT_COUNTS'] * np.sqrt(hdu[b1].data['OPT_COUNTS_IVAR']))
        aa=hdu[b1].data['OPT_WAVE'] * hdu[r1].data['OPT_WAVE']        
    except:
        rSN=0
        bSN=0

    return rSN, bSN                


#######################################################
# GRAB ARC FITS RESIDUALS
#######################################################
def get_arc_fit_residuals(slits,data_dir):

    for det in np.arange(1,9,1):   

        wavefile = data_dir + '/Masters/MasterWaveCalib_A_1_0'+str(det)+'.json'
        jdict = linetools.utils.loadjson(f)
        mdet = slits['det'] == det
       # for i in slits[mdet]
            
        
    return slits               


                                
#names=('bext', 'rext','bdet','rdet', 'bspat','rspat','xpos'))

#######################################################
# MATCH SLITS BASED ON RA
#######################################################
def spec1d_match_red_blue(aslits):

    # ADD DETECTOR NAME
    det = []
    for obj in aslits:
        tmp =   obj['name'].split('DET')
        tdet = int(tmp[1])
        det.append(tdet)
    aslits['det'] = det


    # ***FOR THE MOMENT, REMOVE SERENDIPS
    m=aslits['name'] == 'SERENDIP'
    slits = aslits[~m]
    
    # MATCH RED TO BLUE VIA RA/DEC
    mb = slits['det'] <=4
    mr = slits['det'] >4
    rslits = slits[mr]
    bslits = slits[mb]

    n=0

    # SEARCH ON BLUE FIRST
    for obj in bslits:

        mtc = (obj['objra'] == rslits['objra'])
        if (np.sum(mtc)==1):
            robj = rslits[mtc]
            if (obj['objdec'] != robj['objdec']):
                print('DEC does not match RA!')

            # START ARRAY
            if (n==0):
                matches = Table([[obj['name']],[robj['name']],[obj['det']],[robj['det']],\
                            [obj['objra']],[obj['objdec']],[obj['objname']],[obj['maskdef_id']],[obj['slit']]], \
                            names=('bname', 'rname','bdet','rdet', 'objra','objdec','objname','maskdef_id','xpos'))
            if (n > 0):
                matches.add_row((obj['name'],robj['name'],obj['det'],robj['det'],\
                                 obj['objra'],obj['objdec'],obj['objname'],obj['maskdef_id'],obj['slit']))
            n=n+1

        # NO RED MATCH
        if (np.sum(mtc)==-11): 
        #if (np.sum(mtc)==0):        

            if (n==0):
                matches = Table([[obj['name']],['-1'],[obj['det']],[-1],\
                             [obj['objra']],[obj['objdec']],[obj['objname']],[obj['maskdef_id']],[obj['slit']]], \
                             names=('bname', 'rname','bdet','rdet', 'objra','objdec','objname','maskdef_id','xpos'))
            if (n > 0):
                matches.add_row((obj['name'],'-1',obj['det'],-1,\
                                 obj['objra'],obj['objdec'],obj['objname'],obj['maskdef_id'],obj['slit']))
            n=n+1


    # SEARCH RED OBJECTS FOR NON-MATCHES IN BLUE
    for obj in rslits:

        mtc = (obj['objra'] == bslits['objra'])
        #if (np.sum(mtc)==0):

         #   matches.add_row(('-1',obj['name'],-1,obj['det'],\
         #                        obj['objra'],obj['objdec'],obj['objname'],obj['maskdef_id'],obj['slit']))
         #   n=n+1

    return matches,n


#####################################################
# CREATE SLIT TABLE WHICH INCLUDES FITS FOR FLEXTURE
# CORRECTION, ARC/SKY WAVE RESID 
######################################################
def initialize_slits(nslits, rbext,hdu):

    slits  = mk_slit_table(nslits, hdu)

    for i in np.arange(0,nslits,1):
        
        
        slits['rname'][i] = rbext['rname'][i][0]
        slits['bname'][i] = rbext['bname'][i]
        slits['rdet'][i] = rbext['rdet'][i]
        slits['bdet'][i] = rbext['bdet'][i]

        slits['xpos'][i] = rbext['xpos'][i]
        slits['objra'][i] = rbext['objra'][i]
        slits['objdec'][i] = rbext['objdec'][i]

        slits['slitname'][i] = rbext['objname'][i]
        slits['maskdef_id'][i] = rbext['maskdef_id'][i]
        
        r =slits['bname'][i] 
        hdr = hdu[r].header
#        slits['rms_arc_r'][i] = hdr['WAVE_RMS']

        rSN, bSN = calc_rb_SN(slits['rname'][i],slits['bname'][i], hdu)
        slits['rSN'][i] = rSN
        slits['bSN'][i] = bSN

        b=slits['bname'][i]
        try:
            slits['ypos'][i]      = np.min(hdu[b].data['OPT_WAVE'])
            hd = hdu[b].header
#            slits['rms_arc_b'] = hd['WAVE_RMS']
        except:
            slits['ypos'][i] = -1.
            slits['rSN'][i] = -1.
            slits['bSN'][i] = -1.

    return slits



#####################################################
# CREATE SLIT TABLE WHICH INCLUDES FITS FOR FLEXTURE
# CORRECTION, ARC/SKY WAVE RESID 
######################################################
def create_slit_table(hdu,data_dir,txt):

    header = hdu[0].header   
    
    # READ PYPEIT TXT FILE AND ADD DET COLUMN
    spec1d_txt =  ascii.read(txt,delimiter='|',data_start=1,header_start=0)
    spec1d_txt.rename_column('col0', 'det')

    # MATCH SLITS
    rb_ext, nslits = spec1d_match_red_blue(spec1d_txt)

    
    slits = initialize_slits(nslits,rb_ext,hdu)

    #slits = get_arc_residuals(slits,data_dir)

    return slits,nslits




