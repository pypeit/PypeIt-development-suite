;-----------------------------------------------------------------------
pro deimos_throughput_web_pages, CLOBBER=clobber, TEST=test
;-----------------------------------------------------------------------
;+
; NAME:
;	DEIMOS_THROUGHPUT_WEB_PAGES
;
; PURPOSE:
;	This procedure will generate web pages for DEIMOS throughput
;
; CATEGORY:
;	Instrument performance monitoring
;
; CALLING SEQUENCE:
;       deimos_throughput_web_pages
;
; KEYWORDS:
;       CLOBBER:        if set, overwrite previously-generated plots;
;       otherwise, skip over them.
;
;       TEST:   run in test mode; only generate data for files in
;       directory "test".
;
; RESTRICTIONS:
;       Must run as dmoseng
;
; INPUT FILES:
;       $DEIMOS_CALIBS/<grating>/sens*.fits*
;
; OUTPUT FILES:
;       $DEIMOS_CALIBS/<grating>/doc/EFF_sens_<grating>_DDMMMYYY_NNN.pdf
;       $DEIMOS_CALIBS/<grating>/doc/EFF_sens_<grating>_DDMMMYYY_NNN.png
;       $DEIMOS_CALIBS/<grating>/doc/EFF_sens_<grating>_DDMMMYYY_NNN.txt
;       $DEIMOS_CALIBS/<grating>/doc/ZPANG_sens_<grating>_DDMMMYYY_NNN.pdf
;       $DEIMOS_CALIBS/<grating>/doc/ZPANG_sens_<grating>_DDMMMYYY_NNN.png
;       $DEIMOS_CALIBS/<grating>/doc/ZPANG_sens_<grating>_DDMMMYYY_NNN.txt
;       $DEIMOS_CALIBS/<grating>/doc/ZPPIX_sens_<grating>_DDMMMYYY_NNN.pdf
;       $DEIMOS_CALIBS/<grating>/doc/ZPPIX_sens_<grating>_DDMMMYYY_NNN.png
;       $DEIMOS_CALIBS/<grating>/doc/ZPPIX_sens_<grating>_DDMMMYYY_NNN.txt
;
; PROCEDURE:
;
; REQUIREMENTS:
;       IDL must be invoked using the command
;               ~dmoseng/bin/do_deimos_throughput
;       in order to properly configure the IDL PATH
;
; EXAMPLE:
;       1) Update DEIMOS throughput web pages:
;               deimos_throughput_web_pages
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2009-Oct-23	GDW	Original version
;-
;-----------------------------------------------------------------------

;; specify set of wavelengths at which to measure change over time,
;; plus an array of widths and an array to hold measurements...
lambda_eff = [5000., 6200., 8000., 9000.]
dlambda_eff = [500., 500., 500., 500.]
efficiency = fltarr(n_elements(lambda_eff))

;; build list of directories...
caldir = getenv('DEIMOS_CALIBS')
if caldir eq '' then message, 'DEIMOS_CALIBS envar not defined'
subdirs = ['1200G','600ZD','830G','900ZD']

;; check for test mode...
if keyword_set(TEST) then subdirs = ['test']

;; define contents of a structure
foo = {PARAMS, infile:'', $
       dataset:'', $
       detail:'', $
       fig_eff:'', $
       fig_zp_pix:'', $
       fig_zp_ang:'', $
       fig_eff_pdf:'', $
       fig_zp_pix_pdf:'', $
       fig_zp_ang_pdf:'', $
       tab_eff:'', $
       tab_zp_pix:'', $
       tab_zp_ang:'', $
       date:'', $
       jd:0.d0, $
       std_name:'', $
       ra:'', $
       dec:'', $
       airmass:0., $
       blocking:'', $
       spec_bin:'', $
       grating:'', $
       cenlam:'', $
       slit_width:0., $
       conditions:'', $
       lambda_eff:lambda_eff, $
       dlambda_eff:dlambda_eff, $
       efficiency:efficiency }

;; define structure for summary...
summary = { efficiency_current_plot:'', $
            efficiency_current_pdf:'', $
            efficiency_current_tab:'', $
            efficiency_current_bintab:'', $
            eff_vs_wavelength_plot:'', $
            eff_vs_wavelength_pdf:'', $
            eff_vs_wavelength_tab:'', $
            eff_vs_wavelength_bintab:'', $
            eff_vs_time_plot:'', $
            eff_vs_time_pdf:'', $
            eff_vs_time_tab:'', $
            eff_vs_time_bintab:''}

;; loop over directories (gratings)...
ndirs = n_elements(subdirs)
for i=0,ndirs-1 do begin

    grating = subdirs[i]
    message, 'processing subdirectory '+grating, /info

    ;; get list of files in directory...
    subdir = getenv('DEIMOS_CALIBS')+'/'+grating
    files = file_search( subdir, 'sens*.fits*', count=nfiles)

    ;; no files ==> go to next subdir...
    if nfiles eq 0 then begin
        message, 'no files in subdir '+subdir, /info
        continue
    endif 

    ;; create output directory as needed...
    outdir = subdir + '/doc'
    if ~ file_test( outdir, /dir) then file_mkdir, outdir

    ;; Sort by date
    all_dates = dblarr(nfiles)
    all_grating = strarr(nfiles)
    for jj=0L,nfiles-1 do begin
        meta = xmrdfits(files[jj],1, /silent)
        ;; Convert to Julian
        all_dates[jj] = x_setjdate(meta.date)
        all_grating[jj] = strtrim(meta.grating,2)

    endfor

    ;; verify that all files in here belong to just one grating...
    uni = uniq(all_grating)
    if n_elements(uni) GT 1 then begin
        message, 'directory '+subdir+' contains files from different gratings!', /inf
        continue
    endif

    ;; sort files into order from newest to oldest...
    order = reverse(sort(all_dates))
    files = files[order]

    ;; generate structure to hold data...
    params = replicate({PARAMS}, nfiles)

    ;; intialize params...
    params.infile = files
    for jj=0L,nfiles-1 do begin
        params[jj].lambda_eff  = lambda_eff
        params[jj].dlambda_eff = dlambda_eff
        params[jj].jd = all_dates[jj]
    endfor

    ;; generate plots...
    deimos_throughput_grating_plots, params, fig_time, OUTDIR=outdir, $
      CLOBBER=clobber

    ;; generate detail pages...
    extn = '.html'
    for jj=0L,nfiles-1 do begin

        ;; grab dataset name for input file...
        istrt = strpos( params[jj].infile, '/', /reverse_search) > 0L
        iend = strpos( params[jj].infile, '.fits')
        dataset = strmid( params[jj].infile, istrt+1, iend-istrt-1)

        ;; build output file name...
        outfile = outdir + '/' + dataset + extn

        params[jj].detail = outfile
        params[jj].dataset = dataset
        if keyword_set(VERBOSE) then begin
            message, 'creating detail page '+outfile, /in
        endif
        deimos_throughput_grating_detail_web_page, params[jj]
    endfor

    ;; generate plots for grating summary...
    deimos_throughput_grating_summary_plots, params, summary, $
      OUTDIR=outdir, $
      CLOBBER=clobber

    ;; generate grating summary page...
    extn = '.html'
    outfile = outdir + '/index.html'
    if keyword_set(VERBOSE) then begin
        message, 'creating summary page '+outfile, /in
    endif 
    deimos_throughput_grating_summary_web_page, outfile, grating, $
      params, summary

;    ;; generate web page with plots...
;    outfile = outdir + '/plots.html'
;    deimos_throughput_grating_web_page, $
;      outfile=outfile, $
;      TITLE='DEIMOS '+grating+' Throughput Measurements', $
;      INDEX='throughput: measurements for ' + grating, $
;      GRATING=grating, $
;      FIG_EFF=params[nfiles-1].fig_eff, $
;      FIG_ZP_ANG=params[nfiles-1].fig_zp_ang
;;;      FIG_TIME=fig_time

endfor 

;; deimos_throughput_master_plot
outfile = getenv('DEIMOS_CALIBS')+'/index.html'
gratings = subdirs
href = subdirs + '/doc/index.html'
deimos_throughput_master_web_page, outfile, gratings, href

end
