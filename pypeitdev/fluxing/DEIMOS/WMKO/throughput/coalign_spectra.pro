;-----------------------------------------------------------------------
function coalign_spectra, wave1, flux1, wave2, flux2, $
                     SHIFT=shift, DISPLAY=display
;-----------------------------------------------------------------------
;+
; NAME:
;	COALIGN_SPECTRA
;
; PURPOSE:
;       Given two similar spectra which may have a wavelength shift,
;       determine the shift and, optionally, shift the second spectrum
;       to align withthe first.
;
; CATEGORY:
;	IPM
;
; CALLING SEQUENCE:
;	coalign_spectra, wave1, flux1, wave2, flux2
;
; INPUTS:
;	wave1:  wavelength array for first spectrum
;	flux1:  flux array for first spectrum
;	wave2:  wavelength array for second spectrum
;	flux2:  flux array for second spectrum
;
; KEYWORD PARAMETERS:
;	SHIFT:	Set this keyword to modify the wavelengths of the
;               second spectrum to align with the first.
;
;       DISPLAY: Set this keyword to generate diagnostic plots
;
; OUTPUTS:
;       This function will return the measured wavelength shift in ;
;       wavelength units.  When ; this value is added to wave2, the
;       spectra will align.
;
; COMMON BLOCKS:
;	BLOCK1:	Describe any common blocks here. If there are no COMMON
;		blocks, just delete this entry.
;
; SIDE EFFECTS:
;	If the SHIFT keyword is set, then the wavlength scale on the
;	second spectrum is modified.
;
; PROCEDURE:
;	You can describe the foobar superfloatation method being used here.
;	You might not need this section for your routine.
;
; EXAMPLE:
;       1) Determine shift between two spectra (returned as delta):
;		delta = coalign_spectra( wave1, flux1, wave2, flux2)
;
;       2) Similar to (1) but shift wavelengths on spectrum2 to align:
;		delta = coalign_spectra( wave1, flux1, wave2, flux2, /shift)
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;	Luca Rizzi, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2012-Feb-07	GDW/LR	Original version
;-
;-----------------------------------------------------------------------

;; mask low-illumination regions...
threshold = 0.01                   ; 1% cutoff
spectral_mask, wave1, flux1, wave1_good, flux1_good, threshold=threshold
spectral_mask, wave2, flux2, wave2_good, flux2_good, threshold=threshold

;; normalize the spectra by subtracting a B-spline fit..
flux1_norm = gdw_normalize( wave1_good, flux1_good)
flux2_norm = gdw_normalize( wave2_good, flux2_good)

;; determine the overlapping wavelength region...
wmin = min( [wave1_good, wave2_good], max=wmax)

;; determine median pixel separation...
n = n_elements(wave1_good)
buf = wave1_good[1:n-1] - wave1_good[0:n-2]
dx = median(buf)

;; set up a new wavelength scale from wmin to wmax with pixel size
;; dx...
npix = nint((wmax-wmin)/dx)
dx = (wmax-wmin)/(npix-1)
wave_common = findgen(npix)*dx + wmin

;; re-sample the spectra on the new grid...
flux1_common = interpol( flux1_norm, wave1_good, wave_common)
flux2_common = interpol( flux2_norm, wave2_good, wave_common)

;; determine the lag...
n = n_elements(wave_common)
nlag = 700
lag = indgen(nlag)
nlag = n_elements(lag)
lag = lag - lag[nlag/2]
result = C_CORRELATE(flux1_common, flux2_common, lag)

;; measure the peak (simple way)...
top = max( result, m)
mshift = lag[m]

;;; measure the peak (better way)...
;cut = mean(result)
;result2 = result - cut
;bad = where( result2 lt 0.)
;result2[bad] = 0.
;sum1 = total( result2 * findgen(nlag))
;sum2 = total( findgen(nlag))
;m2 = nint(sum1/sum2)
;mshift2 = lag[m2]
;help, m, mshift
;help, m2, mshift2

;; compute the shift...
if mshift gt 0 then begin
    delta = wave_common[0] - wave_common[mshift]
endif else begin
    delta = wave_common[-mshift] - wave_common[0]
endelse 

;; modify wave2 if requested...
if keyword_set(SHIFT) then wave2 += delta

;; optional plots...
if keyword_set(DISPLAY) then begin

    ;; define useful colors...
    red = gdwcolor('red')
    white = gdwcolor('white')
    green = gdwcolor('green')

    ;; plot normalized and resampled spectra...
    wset, 1
    plot, wave_common, flux1_common, color=white, title='Normalized spectra (white=template, red=target)'
    oplot, wave_common, flux2_common, color=red

    ;; plot the lag...
    wset, 2
    plot, lag, result, xtitle='lag', ytitle='cross-correlation strength', title='delta='+stringify(delta)
    oplot, [mshift,mshift], !y.crange, linestyle=1
    ;; oplot, lag, result2+cut, linestyle=1
    ;; oplot, [mshift2,mshift2], !y.crange, linestyle=1
endif 

return, delta

end
