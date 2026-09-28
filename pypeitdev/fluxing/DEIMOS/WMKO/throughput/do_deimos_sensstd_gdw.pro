;------------------------------------------------------------------------
pro do_deimos_sensstd_gdw, INPUT=input, CLOBBER=clobber, $
  TARGNAME=targname, GRATENAM=gratenam
;------------------------------------------------------------------------
;+
; NAME:
;	DO_DEIMOS_SENSSTD_GDW
;
; PURPOSE:
;	This routine will take output from deimos_throughput.pro 
;       and convert it to a throughput measurement using the Xavier
;       procedure. 
;
; CATEGORY:
;	Instrumentation
;
; CALLING SEQUENCE:
;	DO_DEIMOS_SENSTD_GDW, INPUT=input, OUTPUT=output
;
; INPUT KEYWORDS:
;
;       INPUT: list of DEIMOS throughput images to extract.  Default
;       is to scan all of the directories in the extract/ partition
;       and process all data.
;
;       TARGNAME: set this keyword to limit reduction to the specified
;       target name.
;
;       GRATENAM: set this keyword to limit reduction to the specified
;       grating name.
;
;       CLOBBER: if set, then overwrite existing files (default: no
;       overwrite) 
;
;       DEBUG: set this keyword to generate additional output for
;       debugging purposes.
;
; REQUIREMENTS:
;       IDL must be invoked using the command
;               ~dmoseng/bin/do_deimos_throughput
;       in order to properly configure the IDL PATH
;
; INPUT FILES:
;       $DEIMOS_THRU_DIR/extract/<grating>/YYYYmmmDD_dNNNN_NNNN.sav
;               FITS binary table with two extensions:
;                       - METADATA has keyword/value pairs
;                       - SPECTRUM has wavelength and flux
; OUTPUT FILES:
;       $DEIMOS_THRU_DIR/calib
;       $DEIMOS_THRU_DIR/calib/log/do_deimos_sensstd_gdw.YYYY-MMM-DD-HH:MM:SS.log
;               logfile from execution of this routine
;               
;
; PROCEDURE:
;
; INVOKED SUBROUTINES:
;       deimos_throughput
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2014-Feb-18	GDW	Original version
;-
;------------------------------------------------------------------------

;; set up the output directory based on envar...
thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq '' then begin
    message, 'DEIMOS_THRU_DIR is not defined -- abort!'
endif
outdir = thru_dir + '/calibs'

;; open logfile...
logdir = outdir + '/log'
file_mkdir, logdir
now = gdwtimestamp()
logfile = logdir + '/do_deimos_sensstd_gdw.' + now + '.log'
openw, logunit, logfile, /get_lun
plog, logunit, 'Started'

;; build input and output file lists...
if keyword_set(INPUT) then begin

    infiles = input
    explicit_input_files = 1B

endif else begin
    
    ;; if no input, we scan the directories for all available FITS files...
    dir = getenv('DEIMOS_THRU_DIR') + '/extract'
    infiles = file_search(dir,'*.sav')
    infiles2 = file_search(dir,'*.sav.gz')
    if n_elements(infiles2) gt 0 then infiles = [infiles,infiles2]
    explicit_input_files = 0B

endelse
n_in = n_elements(infiles)

;; loop over input files...
for i=0,n_in-1 do begin

    in = infiles[i]
    plog, logunit, 'Processing file '+in

    ;; verify existence of input file...
    if file_test( in, /regular) then begin
        if keyword_set(DEBUG) then begin
            plog, logunit, 'verified input file '+in+' exists'
        endif 
    endif else begin
        errmsg = 'specified input file '+in+' not found; skipping'
        plog, logunit, in+' '+errmsg
        continue
    endelse

    ;; if .gz file, then verify that it's OK...
    if keyword_set(GZIP_TEST) and stregex( in, '.gz$', /bool) then begin
        plog, logunit, 'Testing GZIP format file'
        cmd = 'gzip -t '+in
        spawn, cmd, buf, exit=status

        ;; if file is bad then rename it...
        if status ne 0 then begin
            badname = in+'.BAD'
            cmd = '/bin/mv '+in+' '+badname
            spawn, cmd, buf, exit=status
            if status eq 0 then begin
                plog, logunit, 'WARNING: renamed bad GZIP file to '+badname
            endif else begin
                plog, logunit, 'WARNING: unable to rename bad GZIP file to '+badname
            endelse 
            continue
        endif 
    endif

    ;; grab grating and target names...
    extract = read_extract_sav_file( in)
    m = extract.meta
    gratenam2 = m.GRATING
    targname2 = m.STD_NAME

    ;; optional match for targname...
    if keyword_set(TARGNAME) then begin
        if strtrim(TARGNAME,2) ne targname2 then begin
            plog, logunit, 'TARGNAME ('+targname2+') is not '+targname+' -- skipping'
            continue
        endif 
    endif

    ;; optional match for graname...
    if keyword_set(GRATENAM) then begin
        if strtrim(GRATENAM,2) ne gratenam2 then begin
            plog, logunit, 'GRATENAM ('+gratenam2+') is not '+gratenam+' -- skipping'
            continue
        endif 
    endif

    ;; convert to sensitivity standard...
    deimos_sensstd_gdw, in, CLOBBER=clobber, LOGUNIT=logunit

endfor 

plog, logunit, 'Results logged to file '+logfile
plog, logunit, 'Finished'
free_lun, logunit
end
