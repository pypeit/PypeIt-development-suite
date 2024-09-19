import os
from astropy.io import fits



#quasar = 'J0411'
#quasar = 'J0020'
quasar = 'J1342'

if quasar == 'J0411':
    #rawpath = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1222/level_12/01222/'
    boguspath = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1222/level_12/01222_bogus_F100LP/'

    exp_list = []
    for detname in ['nrs1', 'nrs2']:
        # NIRSPEC 3-point dither dither center
        scifile1 = os.path.join(boguspath, 'bogus_jw01222002001_03102_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(boguspath, 'bogus_jw01222002001_03102_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(boguspath, 'bogus_jw01222002001_03102_00003_' + detname + '_rate.fits')
        scifile4 = os.path.join(boguspath, 'bogus_jw01222002001_03104_00001_' + detname + '_rate.fits')
        scifile5 = os.path.join(boguspath, 'bogus_jw01222002001_03104_00002_' + detname + '_rate.fits')
        scifile6 = os.path.join(boguspath, 'bogus_jw01222002001_03104_00003_' + detname + '_rate.fits')
        exp_list += [scifile1, scifile2, scifile3, scifile4, scifile5, scifile6]

    for file in exp_list:
        fits.setval(file, 'FILTER', value='F100LP')

elif quasar == 'J0020':
    boguspath = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0020/calwebb/output/'
    exp_list = []
    for detname in ['nrs1', 'nrs2']:
        # NIRSPEC 3-point dither dither center
        scifile1 = os.path.join(boguspath, 'jw01222012001_03102_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(boguspath, 'jw01222012001_03102_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(boguspath, 'jw01222012001_03102_00003_' + detname + '_rate.fits')
        scifile4 = os.path.join(boguspath, 'jw01222012001_03104_00001_' + detname + '_rate.fits')
        scifile5 = os.path.join(boguspath, 'jw01222012001_03104_00002_' + detname + '_rate.fits')
        scifile6 = os.path.join(boguspath, 'jw01222012001_03104_00003_' + detname + '_rate.fits')

        exp_list += [scifile1, scifile2, scifile3, scifile4, scifile5, scifile6]

        # Make copies of the level1 rate files
        for file in exp_list:
            newfile = os.path.join(boguspath, os.path.basename(file).replace('jw', 'bogus_jw'))
            os.system('cp ' + file + ' ' + newfile)
            fits.setval(newfile, 'FILTER', value='F100LP')

elif quasar == 'J1342':
    boguspath = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_1219/calwebb/output/'
    exp_list = []
    slits = ['S200A1', 'S200A2']
    for disperser in ['140H', '235H']:
        for slit in slits:  
            for detname in ['nrs1', 'nrs2']:
                exp_list = []
                if disperser == '140H':
                    if slit == 'S200A1':
                        scifile1 = os.path.join(boguspath, 'jw01219006001_04101_00001_' + detname + '_rate.fits')
                        scifile2 = os.path.join(boguspath, 'jw01219006001_04101_00002_' + detname + '_rate.fits')
                        scifile3 = os.path.join(boguspath, 'jw01219006001_04101_00003_' + detname + '_rate.fits')
                    elif slit == 'S200A2':
                        scifile1 = os.path.join(boguspath, 'jw01219006001_10101_00001_' + detname + '_rate.fits')
                        scifile2 = os.path.join(boguspath, 'jw01219006001_10101_00002_' + detname + '_rate.fits')
                        scifile3 = os.path.join(boguspath, 'jw01219006001_10101_00003_' + detname + '_rate.fits')
                elif disperser == '235H':
                    if slit == 'S200A1':
                        scifile1 = os.path.join(boguspath, 'jw01219006001_06101_00001_' + detname + '_rate.fits')
                        scifile2 = os.path.join(boguspath, 'jw01219006001_06101_00002_' + detname + '_rate.fits')
                        scifile3 = os.path.join(boguspath, 'jw01219006001_06101_00003_' + detname + '_rate.fits')
                    elif slit == 'S200A2':
                        scifile1 = os.path.join(boguspath, 'jw01219006001_08101_00001_' + detname + '_rate.fits')
                        scifile2 = os.path.join(boguspath, 'jw01219006001_08101_00002_' + detname + '_rate.fits')
                        scifile3 = os.path.join(boguspath, 'jw01219006001_08101_00003_' + detname + '_rate.fits')
            
                exp_list += [scifile1, scifile2, scifile3]

                print('---------------------------------------------------')
                print('Creating bogus files for disperser: ' + disperser + ' and slit: ' + slit + ' and detector: ' + detname)
                print(exp_list)
                print('---------------------------------------------------')
                # Make copies of the level1 rate files
                for file in exp_list:
                    # Create bogus fixed slit observations for all setups
                    newfile_fs = os.path.join(boguspath, os.path.basename(file).replace('jw', 'bogus_FS_jw'))
                    print('cp -f ' + file + ' ' + newfile_fs)
                    print("fits.setval(" + newfile_fs + ", 'EXP_TYPE', value='NRS_FIXEDSLIT')")
                    print("fits.setval(" + newfile_fs + ", 'FXD_SLIT', value=" + slit + ")")
                    print("fits.setval(" + newfile_fs + ", 'SRCTYAPT', value='POINT')")
                    os.system('cp -f ' + file + ' ' + newfile_fs)
                    fits.setval(newfile_fs, 'EXP_TYPE', value='NRS_FIXEDSLIT')
                    fits.setval(newfile_fs, 'FXD_SLIT', value=slit)
                    fits.steval(newfile_fs, 'SRCTYAPT', value='POINT')
                    # Also create bogus fixed slit observations with F100LP filter for 140H disperser
                    if disperser == '140H':
                        newfile_fs_filter = newfile_fs.replace('FS', 'FS_F100LP')
                        print('cp -f ' + newfile_fs + ' ' + newfile_fs_filter)
                        print('fits.setval(' + newfile_fs_filter + ', \'FILTER\', value=\'F100LP\')')
                        os.system('cp -f ' + newfile_fs + ' ' + newfile_fs_filter)
                        fits.setval(newfile_fs_filter, 'FILTER', value='F100LP')


