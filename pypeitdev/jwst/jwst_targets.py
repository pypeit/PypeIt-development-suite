import os
from IPython import embed 

def jwst_spec1d_files(progid, disperser, target, slit):
    """
    Routine to return a list of spec1d filenames for JWST NIRSpec exposures for a given target and disperser and possibly slit.
    
    Parameters
    ----------
    progid : str
        Program ID for the observations. 
    disperser : str
        Name of the disperser.
    target : str
        Name of the target.
    slit : str
        Slit requested, optional default is None.
    
    Return 
    ------
    spec1d_list : list
        List of spec1d filenames.
    """

    uncal_list, redux_path, rawpath_level2 = jwst_targets(progid, disperser, target, slit=slit)
    spec1d_pre = ['spec1d_' + os.path.basename(file).split('_nrs1')[0] + f'_slit_{slit}.fits' for file in uncal_list[0]]
    spec1d_filenames = [os.path.join(redux_path, 'pypeit', 'Science', spec) for spec in spec1d_pre]
    
    return redux_path, spec1d_filenames

def jwst_targets(progid, disperser, target, slit=None):
    """
    Routine to return a list of JWST NIRSpec exposures for a given target and disperser and possibly slit.
    
    Parameters
    ----------
    progid : str
        Program ID for the observations. 
    disperser : str
        Name of the disperser.
    target : str
        Name of the target.
    slit : str
        Slit requested, optional default is None.
        
    Returns
    -------
    exp_list : list
        List of lists of uncalibrated files for each exposure.
    redux_dir : str
        Path to the directory where the data will be reduced.
    rawpath_level2 : str
        Path to the directory where the raw level 2 output files are.
    
    """

    exp_list = []
    detectors = ['nrs1', 'nrs2']

    # If bkg_redux is False, the code will model the sky and the object profile and perform optimal extraction.
    # If bkg_redux is True, the code will difference image and simply boxcar extract (optimal not implemented yet)
    for detname in detectors:
        if '2073' in progid:
            if 'PRISM' == disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2073/level_12/02073/'
                redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02073_CLEAR_PRISM', target)
                calwebb_dir = os.path.join(redux_dir, 'output')
                pypeit_output_dir = os.path.join(redux_dir, 'pypeit')

                #J0252
                if target == 'J0252':
                    uncalfile1 = os.path.join(rawpath_level2, 'jw02073007001_03101_00001_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw02073007001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw02073007001_03101_00003_' + detname + '_uncal.fits')

                elif target == 'J1007':
                    # J1007
                    # NIRSPEC 3-point dither
                    uncalfile1 = os.path.join(rawpath_level2, 'jw02073008001_03101_00001_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw02073008001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw02073008001_03101_00003_' + detname + '_uncal.fits')

                    #uncalfile1 = os.path.join(rawpath_level2, 'jw02073006001_03101_00001_' + detname + '_uncal.fits')
                    #uncalfile2 = os.path.join(rawpath_level2, 'jw02073006001_03101_00002_' + detname + '_uncal.fits')
                    #uncalfile3 = os.path.join(rawpath_level2, 'jw02073006001_03101_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3]) 

        elif '2756' in progid:
            if 'PRISM' in disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2756/level_12/02756/'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02756_CLEAR_PRISM/calwebb'
                output_dir = os.path.join(redux_dir, 'output')
                pypeit_output_dir = os.path.join(redux_dir, 'pypeit')

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw02756001001_03101_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02756001001_03101_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02756001001_03101_00003_' + detname + '_uncal.fits')

                # dither offset
                # uncalfile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_uncal.fits')
                # bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_uncal.fits')
                # bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        if '1133' in progid:
            if 'PRISM' in disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/pypeit'

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_uncal.fits')

                # dither offset
                # uncalfile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_uncal.fits')
                # bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_uncal.fits')
                # bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        elif '2072' in progid:
            if 'PRISM' in disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/2072/level_12'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/2072/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/02027_PRISM/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/02027_PRISM/calwebb/pypeit'

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw02072002001_05101_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02072002001_05101_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02072002001_05101_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        elif '1219' in progid:
            if 'J1342+0928' in target:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_1219/1219/level_12/01219/'
                redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/1219/J1342+0928/')

                if disperser == '140H':
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_04101_00001_' + detname + '_uncal.fits')  # msa_metadata_id  = 1
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_04101_00002_' + detname + '_uncal.fits')  
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_04101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_10101_00001_' + detname + '_uncal.fits') # msa_metadata_id  = 29
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_10101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_10101_00003_' + detname + '_uncal.fits')
                elif disperser == '235H':
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_06101_00001_' + detname + '_uncal.fits') # msa_metadata_id  = 15
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_06101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_06101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_08101_00001_' + detname + '_uncal.fits') # msa_metadata_id  = 16
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_08101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_08101_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])           

        elif '1222' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1222/level_12/01222/'
            if 'J0411-0907' in target:
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1222/J0411-0907/'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03108_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03106_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03104_00003_' + detname + '_uncal.fits')

            elif 'J0020-3653' in target:
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1222/J0020-3653/'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03108_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03106_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03104_00003_' + detname + '_uncal.fits')
                        
            elif 'J1120+0641' in target:
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1222/J1120+0641/'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03108_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03106_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03104_00003_' + detname + '_uncal.fits')                        

            exp_list.append([uncalfile1, uncalfile2, uncalfile3])

        if '1764' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1764/level_12/01764/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/', target)
            if 'J0313-1806' in target:
                #redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J0313-1806/calwebb'
                #output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J0313-1806/calwebb/output'
                #pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J0313-1806/calwebb/pypeit'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03104_00003_' + detname + '_uncal.fits')
                elif disperser == '395H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03106_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03108_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00003_' + detname + '_uncal.fits')
                
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

            elif 'J1007+2115' in target:
                #redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J1007+2115/calwebb'
                #output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J1007+2115/calwebb/output'
                #pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J1007+2115/calwebb/pypeit'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04104_00003_' + detname + '_uncal.fits')
                elif disperser == '395H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04106_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04108_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        if '3526' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/3526/level_12/03526/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/3526/', target)
            #output_dir = os.path.join(redux_dir, 'calwebb', 'output')        
            #pypeit_output_dir = os.path.join(redux_dir, 'calwebb', 'pypeit')
            # All of these are with slit S200A2
            if 'J0410-0139' in target:
                file_list = []
                for ii in range(1,11): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526012001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            elif 'J0038-1527' in target:
                file_list = []
                for ii in range(1,4): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526007001_05101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)  
            elif 'J0252-0503' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526008001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                 
            elif 'J0038-0653' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526001001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            # All of these are with slit S200A1
            elif 'J0313-1806' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526002001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                       
            elif 'J1007+2115' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526005001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                

        if '1117' in progid:
            if 'PRISM' in disperser:
                # PRISM data
                rawpath_level2 = '//Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/level_12/01117'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/pypeit'

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw01117007001_03101_00002_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw01117007001_03101_00003_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw01117007001_03101_00004_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        if '2736' in progid:
            if 'G395M' in disperser:
                # Use islit = 37 for nrs1
                # G395M data
                rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/Raw'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/pypeit'

                # NIRSPEC 3-point dither
                uncalfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_uncal.fits')
                
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

            elif 'G235M' in disperser:
                # Use islit = 38 for nrs1
                # G235M data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/pypeit/'

                # NIRSPEC 3-point dither
                uncalfile1 = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_uncal.fits')
            
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

        if '1671' in progid:
            if 'G395M' in disperser:
                # Use islit = 37 for nrs1
                # G395M data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/Maseda/'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb/pypeit'

                # NIRSPEC 3-point dither
                uncalfile1 = os.path.join(rawpath_level2, 'jw01671001001_03101_00002_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw01671001001_03101_00003_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw01671001001_03101_00004_' + detname + '_uncal.fits')
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])


    return exp_list, redux_dir, rawpath_level2
    
    