
type = 'ESI'
;;
CASE type OF
   'ESI': BEGIN
      flatfile = $
         '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/FlatECH10_1x1_D.fits.gz'
      piximgfile = $
         '/Users/joe/REDUX/esi_redux/Mar_2008/Final/f_esi1044.fits.gz'
      waveimgfile = $
         '/Users/joe/REDUX/esi_redux/Mar_2008/Arcs/ArcECH10_1x1IMG.fits.gz'
      sedg_file = $
         '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/SEdgECH10_1x1.fits.gz'
      tset_slits = mrdfits(sedg_file, 1)
      flat = mrdfits(flatfile,0)
      rn = 2.5
      flatvar = abs(flat - sqrt(2.0)*rn) + rn^2
      flativar =(flatvar GT 0.0)/(flatvar + (flatvar EQ 0.0))
      piximg = mrdfits(piximgfile,3)
      waveimg = mrdfits(waveimgfile,0)
      ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit)
      order_vec = -(lindgen(nslit) +1) + 16L
      ordermask = esi_ordermask(tset_slits)

   END
   'MAGE': BEGIN
      mage_path = '/Users/joe/REDUX/mage_redux/redux/mar09/slit0.70/'
      orderfile = mage_path + 'Orders.fits'
      tset_slits = mrdfits(orderfile,1)
      ordermask = mrdfits(orderfile,0)
      piximgfile = mage_path + 'piximg_flats.fits'
      piximg = mrdfits(piximgfile)
      flatfile ='/Users/joe/pypeit_vs_lowredux/MAGE/mage1024.fits.gz'
      mage_proc, flatfile, flat, flativar
      ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit)
      order_vec = -(lindgen(nslit) +1) + 21L
      ordermask = mage_ordermask(tset_slits)
      mage_flatfile = mage_path + 'Flat/Pixflat.fits'
      mage_flat = mrdfits(mage_flatfile,0)
   END
   'GNIRS': BEGIN
      gnirs_path = '/Users/joe/GN-LP-7/redux/J1335+3533/'
      superflatfile = $
         '/Users/joe/GN-LP-7/redux/J1335+3533/superflat-N20160127S0413-428.fits'
      order_vec = [3, 4, 5, 6, 7, 8]
      tset_slits =mrdfits(superflatfile,1)
      ordermask = long_slits2mask(tset_slits)
      ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit)
      ordermask[WHERE(ordermask GT 0)] = ordermask[WHERE(ordermask GT 0)] + 2L
      flatfile = gnirs_path + 'gnirs_testflat.fits'
      flat = mrdfits(flatfile)
      rn = 12.0
      flatvar = abs(flat - sqrt(2.0)*rn) + rn^2
      flativar =(flatvar GT 0.0)/(flatvar + (flatvar EQ 0.0))
      piximgfile = gnirs_path + 'gnirs_piximg.fits'
      piximg = mrdfits(piximgfile,0)
      gnirs_flat =mrdfits(superflatfile,0)
      gnirs_proc, '/Users/joe/GN-LP-7/Raw/cN20160127S0409.fits',sciimg_raw
   END
   ELSE: PRINT, 'Unrecognized instrument'
ENDCASE

dims = size(flat, /dimens)
nx = dims[0]
ny = dims[1]

midwidth = tset_slits[1].coeff[0, *] - tset_slits[0].coeff[0, *]


CHK = 1
pixelflat = fltarr(nx, ny) + 1.0D 
illum_flat = fltarr(nx,ny) + 1.0D
model_img = fltarr(nx,ny) + 0.0D
smooth_fit = fltarr(nx,ny)
if (NOT keyword_set(npoly1)) then npoly1 = 3
FOR slit = 0L, nslit-1 DO BEGIN
   igood = where(ordermask EQ order_vec[slit] AND piximg GT 0, ngood)
   
   npercol = floor(float(ngood) / ny) > 1
   npoly = npoly1 < ceil(npercol / 10.)
   npoly = npoly > 1
   
   
  long_flatfield_specillum,piximg[igood], ximg[igood] $
                           , flat[igood], spec_set $
                           , illum_set $
                           , xsamp = 1.5 $
                           , SLITSAMP=slitsamp $
                           , invvar = flativar $ 
                           , slitwidth = midwidth[slit-1] $
                           , modelfit = modelfit $
                           , ybkpt = ybkpt $
                           , npoly = npoly, CHK = CHK $
                           , PIXFIT = 1, smoothfit=smoothfit
  qgood = modelfit GT 0
  pixelflat[igood] = qgood*flat[igood]/(modelfit + (qgood EQ 0))
  illum_flat[igood] = bspline_valu(ximg[igood], illum_set)
  model_img[igood] =  modelfit
  smooth_fit[igood] = smoothfit
ENDFOR

;;
full_flat = illum_flat*pixelflat
sciimg = sciimg_raw
minval = 0.5
divideflat, sciimg, full_flat, minval = minval

END
