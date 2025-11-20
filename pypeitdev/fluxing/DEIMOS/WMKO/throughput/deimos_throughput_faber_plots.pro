pro deimos_throughput_faber_plots, infiles, psfiles, cenlam0

nfiles = n_elements(infiles)
charsize = 1.5
j = 0

;; compile a list of all files matching given central wavelength...
for i=0,nfiles-1 do begin

    ;; define files...
    infile = infiles[i]
    psfile = psfiles[i]

    ;; read data...
    meta = mrdfits( infile, 1, /silent)
    spec = mrdfits( infile, 2, /silent)
    wavelength = spec.wavelength
    flux = spec.counts
    cenlam = nint(meta.central_wave)

    ;; skip data we do not want to compare...
    if cenlam ne cenlam0 then begin
        print, 'Skipping '+infile+' because cenlam='+stringify(cenlam)
        continue
    endif 

    ;; increment counter...
    j += 1

    ;; shift spectrum as needed...
    if j eq 1 then begin

        wavelength0 = wavelength
        flux0 = flux
        infile0 = infile

        ;; fit background with B-spline...
        nord = 4
        fullbkpt = bspline_bkpts( wavelength0, nord=nord, nbkpts=15)
        sset = bspline_iterfit( wavelength0, flux0, nord=nord, $
                                maxiter=0, yfit=fit0, fullbkpt=fullbkpt)

    endif else begin

        ;; sanity check...
        n0 = n_elements(wavelength0)
        n = n_elements(wavelength)
        if n ne n0 then begin
            print, infile, ' has ', n, ' elements'
            print, infile0, ' has ', n0, ' elements'
            message, 'conflict between '+infile+' and '+infile0, /info
            continue
        endif 
        
        coalign_spectra, wavelength0, flux0, wavelength, flux
    endelse 

    
    ;; build title...
    title = meta.std_name + ' ' $
            + meta.grating + '@' $
            + stringify(cenlam) + 'A+' $
            + meta.blocking + ' ' $
            + meta.date

    ;; generate plot...
    set_plot, 'ps'
    device, filename=psfile, /landscape, /color
    plot, wavelength0, alog10(flux0), /nodata, $
          font=1, charsize=charsize, $
          yrange=[2.0,3.5], $
          xtitle='Wavelength [A]', ytitle='log10(Flux)', title=title
    loadct, 13
    oplot, wavelength0, alog10(flux0), color=255
    oplot, wavelength, alog10(flux)

    ;; add filename...
    BREAKNAME, infile, Dirname, Rootname, Extn
    if j eq 1 then begin
        label = rootname
        rootname0 = rootname
    endif else begin
        label = rootname + ' / ' + rootname0
    endelse 
    rlabel, label, charsize=charsize, font=1

    device, /close

endfor 

end
