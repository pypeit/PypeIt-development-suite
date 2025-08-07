import os
from IPython import embed 

def jwst_spec1d_files(progid, disperser, target, slit=None, source=None):
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
    source : str
        Source name, optional default is None.
    
    Return 
    ------
    spec1d_list : list
        List of spec1d filenames.
    """
    if slit is None and source is None:
        raise ValueError("Either 'slit' or 'source' must be specified.")
    elif slit is not None and source is not None:
        raise ValueError("Only one of 'slit' or 'source' can be specified.")


    uncal_list, redux_path, rawpath_level2 = jwst_targets(progid, disperser, target, slit=slit)
    suffix = f'_slit_{slit}.fits' if slit is not None else f'_source_{source}.fits'
    spec1d_pre = ['spec1d_' + os.path.basename(file).split('_nrs1')[0] + suffix for file in uncal_list[0]]
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
        List of lists of uncalibrated files for each exposure. exp_list[0] contains 
        all of the exposure files for nrs1 and exp_list[1] contains all of the exposure files for nrs2.
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
        if '3543' in progid:
            if 'G395M' == disperser:
                ## BHstar object
                
                if target == 'BHstar':
                    rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/3543/level_12/03543/'
                    redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/3543/', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw03543001001_07101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw03543001001_07101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw03543001001_07101_00004_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3]) 
        
        if '2073' in progid:
            if 'PRISM' == disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2073/level_12/02073/'
                #redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02073_CLEAR_PRISM', target)
                redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/2073/', target)

                #J0252
                if target == 'J0252-0503':
                    uncalfile1 = os.path.join(rawpath_level2, 'jw02073007001_03101_00001_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw02073007001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw02073007001_03101_00003_' + detname + '_uncal.fits')

                elif target == 'J1007+2115':
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
                if disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_03101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_03101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_03101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_09101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_09101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_09101_00003_' + detname + '_uncal.fits')
                elif disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_05101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_05101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_05101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_07101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_07101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_07101_00003_' + detname + '_uncal.fits')                        

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

        if '9180' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/9180/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/9180/', target)
            file_list = []
            if 'J2356+0017' in target:
                if disperser == '140H': 
                    prefix = 'jw09180041001_03102_000' if slit == 'S200A1' else 'jw09180041001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180041001_03107_000' if slit == 'S200A1' else 'jw09180041001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J1732+6531' in target:
                if disperser == '140H': 
                    prefix = 'jw09180012001_05101_000' if slit == 'S200A1' else 'jw09180012001_07101_000'
                    indx_range = range(1,6)                    
                elif disperser == '235H':
                    prefix = 'jw09180012001_11101_000' if slit == 'S200A1' else 'jw09180012001_09101_000'
                    indx_range = range(1,6)                    
                elif disperser == '395M':
                    prefix = 'jw09180049001_03102_000'
                    indx_range = range(1,4)
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
            if 'J1429-0104' in target:
                if disperser == '140H': 
                    prefix = 'jw09180031001_03102_000' if slit == 'S200A1' else 'jw09180031001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180031001_03108_000' if slit == 'S200A1' else 'jw09180031001_03106_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,3)           
            if 'J1428+0454' in target:        
                if disperser == '140H': 
                    prefix = 'jw09180029001_03102_000' if slit == 'S200A1' else 'jw09180029001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180029001_03107_000' if slit == 'S200A1' else 'jw09180029001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J1450-0144' in target:        
                if disperser == '140H': 
                    prefix = 'jw09180022001_03102_000' if slit == 'S200A1' else 'jw09180022001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180022001_03107_000' if slit == 'S200A1' else 'jw09180022001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)            
            if 'J1609+5328' in target:    
                if disperser == '140H': 
                    prefix = 'jw09180009001_04101_000' if slit == 'S200A1' else 'jw09180009001_06101_000'
                    indx_range = range(1,4)
                elif disperser == '235H':
                    prefix = 'jw09180009001_10101_000' if slit == 'S200A1' else 'jw09180009001_08101_000'
                    indx_range = range(1,4)                    
                elif disperser == '395M':
                    prefix = 'jw09180050001_04101_000'
                    indx_range = range(1,3)
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
            if 'J1440+0019' in target: 
                # Missing images in mast 
                if disperser == '140H': 
                    prefix = 'jw09180037001_03102_000' if slit == 'S200A1' else 'jw09180037001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180037001_03108_000' if slit == 'S200A1' else 'jw09180037001_03106_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,3)

            for ii in indx_range: 
                file_list.append(os.path.join(rawpath_level2, prefix + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
            exp_list.append(file_list)


        if '1967' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1967/level_12/01967/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/1967/', target)
            # All of these are with slit S200A2
            file_list = []
            if 'J2255+0251' in target:
                for ii in range(1,4): 
                    file_list.append(os.path.join(rawpath_level2, 'jw01967012001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            elif 'J2236+0032' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967011001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            elif 'J0844-0052' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967002001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)            
            elif 'J0844-0132' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967003001_05101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)   
            elif 'J0918+0139' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967005001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list) 
            elif 'J1425-0015' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967008001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                                   
            elif 'J1512+4422' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967009001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      
            elif 'J1525+4303' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967010001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      
            elif 'J1146+0124' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967006001_04102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      
            elif 'J1146-0005' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967007001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)     
            elif 'J0217-0208' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967001001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                     
            elif 'J0911+0152' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967004001_05101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
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
    
    