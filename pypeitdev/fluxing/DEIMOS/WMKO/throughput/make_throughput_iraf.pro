pro make_throughput_iraf, IN=in, OUTFILE=outfile, $
  SIDE=side, GRATENAM=gratenam, CENLAM=cenlam, $
  EXPTIME=exptime, DISPLAY=display

n_truncate = 4

;; count images...
n_in = n_elements(in)

for i=0,n_in-1 do begin

    infile = in[i]
    buf = mrdfits( infile)
    spec = buf[*,0,0]

    ;; flip as needed...
    if side[i] eq 'R' then spec = reverse(temporary(spec))

    ;; count pixels...
    n_pix = n_elements(spec)

    ;; generate wavelengths...
    if side[i] eq 'B' then begin
        wave = deimos_get_wavelengths( gratenam, cenlam, n_pix, /BLUE)
    endif else begin
        wave = deimos_get_wavelengths( gratenam, cenlam, n_pix, /RED)
    endelse

    ;; concatenate spectra...
    if i eq 0 then begin
        wavelength = wave[0:n_pix-n_truncate-1]
        counts = spec[0:n_pix-n_truncate-1]
    endif else begin
        wavelength = [wavelength,wave[n_truncate:n_pix-1]]
        counts = [counts,spec[n_truncate:n_pix-1]]
    endelse

endfor 

if keyword_set(DISPLAY) or keyword_set(PSFILE) then begin
    plot, wavelength, counts, /yno, $
          xtitle='Wavelength (A)', ytitle='Flux [e-/px]', $
          title='Spectrum', font=1
endif

;; divide the counts by the exposure time...
counts = counts/exptime

;; optional output file...
if keyword_set( OUTFILE) then begin
    n = n_elements(wavelength)
    line = {wavelength:0., counts:0.}
    temp = REPLICATE(line,n)
    temp.wavelength = wavelength
    temp.counts     = counts
    sxaddpar, outheader, 'EXTNAME', 'SPECTRUM', 'Name of this extension'
    mwrfits, temp, outfile, outheader
endif 

end
