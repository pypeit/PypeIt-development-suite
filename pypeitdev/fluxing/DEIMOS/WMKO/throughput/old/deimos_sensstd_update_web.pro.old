;-----------------------------------------------------------------------
pro deimos_sensstd_update_web, WVTIME=wvtime
;-----------------------------------------------------------------------

;; process keywords...
if not keyword_set(WVTIME) then wvtime = [5000., 6200., 8000., 9000.]

;; build list of directories...
caldir = getenv('DEIMOS_CALIBS')
subdirs = ['1200G','600ZD','830G','900ZD']

;; loop over directories...
ndirs = n_elements(subdirs)
for i=0,ndirs-1 do begin

    mgrat = subdirs[i]
    message, 'processing subdirectory '+mgrat, /info

    ;; get list of files in directory...
    subdir = getenv('DEIMOS_CALIBS')+'/'+mgrat
    files = file_search( subdir, 'sens_*')
    nfiles = n_elements(files)
    mkhtml_specthru_gdw, files, $
                         TITLE='DEIMOS '+mgrat+' Throughput Measurements', $
                         WVTIME=WVTIME, $
                         OUTPTH=subdir+'/'

endfor 

end
