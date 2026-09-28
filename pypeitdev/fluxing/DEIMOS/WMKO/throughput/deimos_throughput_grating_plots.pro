;-----------------------------------------------------------------------
pro oplot_bar, x, y, dx, dy
;-----------------------------------------------------------------------
thick = 10.0
oplot, [x], [y], psym=7, symsize=2.0, thick=thick
xbar = [-dx,-dx,-dx,dx,dx,dx]
ybar = [dy,-dy,0.,0.,-dy,dy]
oplot, xbar+x, ybar+y, thick=thick
end

;-----------------------------------------------------------------------
pro write_table, outfile, x, y
;-----------------------------------------------------------------------
; write x,y to text file...
;-----------------------------------------------------------------------

openw, outunit, outfile, /get_lun
n = n_elements(x)
for i=0,n-1 do printf, outunit, x[i], y[i]
free_lun, outunit

end

;-----------------------------------------------------------------------
pro deimos_throughput_grating_plots, params, fig_time, $
  OUTDIR=outdir, VERBOSE=verbose, CLOBBER=clobber
;-----------------------------------------------------------------------
; Purpose:
;       generate plots and data tables of zero points and efficiencies
;-----------------------------------------------------------------------

nfiles = n_elements( params)
template = '/tmp/temp.%.ps'
angstrom = STRING(197B)

;; define plot params...
csize = 1.5
xmargin=[3,1]
ymargin=[1,1]

;; loop over files...
for i=0,nfiles-1 do begin

    if keyword_set(VERBOSE) then message, 'processing '+params[i].infile, /info

    ;; Read in the META and THRUPUT data
    meta = xmrdfits( params[i].infile, 1, /silent)
    thru_str = xmrdfits( params[i].infile, 2, /silent)

    ;; grab dataset name for input file...
    istrt = strpos( params[i].infile, '/', /reverse_search) > 0L
    iend = strpos( params[i].infile, '.fits')
    root = strmid( params[i].infile, istrt+1, iend-istrt-1)
    params[i].dataset = root

    ;; generate PNG filenames...
    extn = '.png'
    params[i].fig_zp_pix = OUTDIR+'/ZPPIX_'+root+extn
    params[i].fig_zp_ang = OUTDIR+'/ZPANG_'+root+extn
    params[i].fig_eff   = OUTDIR+'/EFF_'+root+extn

    ;; generate PDF filenames...
    extn = '.pdf'
    params[i].fig_zp_pix_pdf = OUTDIR+'/ZPPIX_'+root+extn
    params[i].fig_zp_ang_pdf = OUTDIR+'/ZPANG_'+root+extn
    params[i].fig_eff_pdf   = OUTDIR+'/EFF_'+root+extn

    ;; generate table names...
    extn = '.txt'
    params[i].tab_zp_pix = OUTDIR+'/ZPPIX_'+root+extn
    params[i].tab_zp_ang = OUTDIR+'/ZPANG_'+root+extn
    params[i].tab_eff   = OUTDIR+'/EFF_'+root+extn

    ;; store useful stuff...
    params[i].date       = meta.date
    params[i].std_name   = meta.std_name
    params[i].ra         = meta.ra
    params[i].dec        = meta.dec
    params[i].airmass    = meta.airmass
    params[i].blocking   = meta.blocking
    params[i].spec_bin   = meta.spec_bin
    params[i].grating    = strtrim(meta.grating,2)
    params[i].cenlam     = stringify(round(meta.central_wave))
    params[i].slit_width = meta.slit_width
    params[i].conditions = meta.conditions

    ;; measure efficiency.  NOTE: can't pass params[i] by reference,
    ;; so to get values back we must copy element [i] into a temporary
    ;; variable...
    buf = params[i]
    deimos_throughput_get_efficiency, thru_str, buf, VERBOSE=verbose
    params[i] = buf

    ;; create a title for all plots...
    title = meta.std_name + ' ' + meta.date + $
            ' X=' + stringify(meta.airmass,'(f4.2)') + $
            ' Bin=' + stringify(meta.spec_bin,'(i2)') + $
            ' Grating=' + meta.grating

    ;;----------------------------------------
    ;; zp pixel plot...
    ;;----------------------------------------
    if (keyword_set(CLOBBER) || $
        ~ file_test(params[i].fig_zp_pix) || $
        ~ file_test(params[i].fig_zp_pix_pdf)) then begin
        psfile = mktemp(template)
        psland, psfile
        DEVICE, SET_FONT='Helvetica', /TT_FONT  
        xrng = [min(thru_str.wav, max=xmax), xmax]
        yrng = [min(thru_str.zp_pix, max=ymax), ymax*1.05]
        plot, thru_str.wav, thru_str.zp_pix, xrange=xrng, $
              yrange=yrng, $
              charsize=csize, $
              xstyle=1, ystyle=1, $
              thick=2, xthick=1, ythick=1, $
              xtitle='Wavelength ['+Angstrom+']', $
              ytitle='Zero Point [AB Mag for 1e-/s/pix]', font=0, title=title, $
              xmargin=xmargin, ymargin=ymargin
        id, font=0
        device, /close
        ps2other, psfile, $
                  png=params[i].fig_zp_pix, $
                  pdf=params[i].fig_zp_pix_pdf, $
                  verbose=verbose, /delete
    endif else begin
        if keyword_set(VERBOSE) then begin
            message, root+' zp_pix plots exist; skipping', /info
        endif 
    endelse 

    ;;----------------------------------------
    ;; zp ang plot...
    ;;----------------------------------------
    if (keyword_set(CLOBBER) || $
        ~ file_test(params[i].fig_zp_ang) || $
        ~ file_test(params[i].fig_zp_ang_pdf)) then begin
        psfile = mktemp(template)
        psland, psfile
        DEVICE, SET_FONT='Helvetica', /TT_FONT  
        DEVICE, /ISOLATIN1
        xrng = [min(thru_str.wav, max=xmax), xmax]
        yrng = [min(thru_str.zp_ang, max=ymax), ymax*1.05]
        plot, thru_str.wav, thru_str.zp_ang, xrange=xrng, yrange=yrng, $
              charsize=csize, $
              xstyle=1, ystyle=1, $
              thick=2, xthick=1, ythick=1, $
              xtitle='Wavelength ['+Angstrom+']', $
              ytitle='Zero Point [AB Mag for 1e-/s/Ang]', font=0, title=title, $
              xmargin=xmargin, ymargin=ymargin
        id, font=0
        device, /close
        ps2other, psfile, $
                  png=params[i].fig_zp_ang, $
                  pdf=params[i].fig_zp_ang_pdf, $
                  verbose=verbose, /delete
    endif else begin
        if keyword_set(VERBOSE) then begin
            message, root+' zp_ang plots exist; skipping', /info
        endif 
    endelse 

    ;;----------------------------------------
    ;; efficiency plot...
    ;;----------------------------------------
    if (keyword_set(CLOBBER) || $
        ~ file_test(params[i].fig_eff) || $
        ~ file_test(params[i].fig_eff_pdf)) then begin
        psfile = mktemp(template)
        psland, psfile
        DEVICE, SET_FONT='Helvetica', /TT_FONT  
        DEVICE, /ISOLATIN1
        xrng = [min(thru_str.wav, max=xmax), xmax]
        yrng = [min(thru_str.eff, max=ymax), ymax*1.05]
        plot, thru_str.wav, thru_str.eff, xrange=xrng, $
              yrange=yrng, $
              charsize=csize, $
              xstyle=1, ystyle=1, $
              thick=2, xthick=1, ythick=1, $
              xtitle='Wavelength ['+Angstrom+']', $
              ytitle='End-to-end Efficiency', $
              title=title, font=0, $
              xmargin=xmargin, ymargin=ymargin
        id, font=0

        ;; overplot efficiency...
        good = where( finite( params[i].efficiency), count)
        if count gt 0 then begin
            for j=0,count-1 do begin
                g = good[j]
                x = params[i].lambda_eff[g]
                y = params[i].efficiency[g]
                dx = params[i].dlambda_eff[g] / 2.
                dy = (yrng[1]-yrng[0])/20.
                oplot_bar, x, y, dx, dy
            endfor
        endif else begin
            if keyword_set(VERBOSE) then begin
                message, 'no good efficiency values', /info
            endif 
        endelse

        device, /close
        ps2other, psfile, $
                  png=params[i].fig_eff, $
                  pdf=params[i].fig_eff_pdf, $
                  verbose=verbose, /delete
    endif else begin
        if keyword_set(VERBOSE) then begin
            message, root+' eff plots exist; skipping', /info
        endif 
    endelse 

    ;;----------------------------------------
    ;; save tables...
    ;;----------------------------------------
    if (keyword_set(CLOBBER) || ~file_test( params[i].tab_zp_pix)) then begin
        write_table, params[i].tab_zp_pix, thru_str.wav, thru_str.zp_pix
    endif else begin
        if keyword_set(VERBOSE) then begin
            message, root+' tab_zp_pix plots exist; skipping', /info
        endif 
    endelse 

    if (keyword_set(CLOBBER) || ~file_test( params[i].tab_zp_ang)) then begin
        write_table, params[i].tab_zp_ang, thru_str.wav, thru_str.zp_ang
    endif else begin
        if keyword_set(VERBOSE) then begin
            message, root+' tab_zp_ang plots exist; skipping', /info
        endif 
    endelse 

    if (keyword_set(CLOBBER) || ~file_test( params[i].tab_eff)) then begin
        write_table, params[i].tab_eff,    thru_str.wav, thru_str.eff
    endif else begin
        if keyword_set(VERBOSE) then begin
            message, root+' tab_eff plots exist; skipping', /info
        endif 
    endelse 

endfor 
end
