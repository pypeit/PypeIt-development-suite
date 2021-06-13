;-----------------------------------------------------------------------
pro deimos_build_throughput_wavelength_tables, infile, outfile
;-----------------------------------------------------------------------

;; infile lists:
;;   1) filename
;;   2) grating
;;   3) cenlam
;;   4) side (blue or red)
readcol, infile, file, gratenam, cenlam, side, format='(a,a,i,a)'

data = { wavelengths: fltarr(4096), $
         gratenam: '', $
         central_wavelength: 0, $
         side: '' $
       }

create = 1B

;; loop over files...
n_files = n_elements(file)
for i=0,n_files-1 do begin

    ;; get wavelengths...
    readcol, file[i], w, format='(f)'

    ;; insert data into struct...
    data.wavelengths        = w
    data.gratenam           = gratenam[i]
    data.central_wavelength = cenlam[i]
    data.side               = side[i]

    ;; write data to file...
    mwrfits, data, outfile, create=create
    create = 0B

endfor 

end
