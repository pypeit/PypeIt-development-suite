pro compare_iraf_files, IN=in, OUT=out, $
  SIDE=side, GRATENAM=gratenam, CENLAM=cenlam

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

endfor 


end
