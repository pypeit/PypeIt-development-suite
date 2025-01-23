;-----------------------------------------------------------------------
pro review_trace, tracefiles, WAIT=wait
;-----------------------------------------------------------------------
;
;-----------------------------------------------------------------------

;; if no files passed, grab all region files in curent directory...
if n_params() eq 0 then tracefiles = file_search( '*.region')

;; count files...
n = n_elements(tracefiles)

;; loop over input files...
for i=0,n-1 do begin

    tracefile = tracefiles[i]

    ;; read header info...
    openr, iunit, tracefile, /get_lun
    line = ''

    readf, iunit, line
    token = STRSPLIT( line, ' ', COUNT=count, /EXTRACT)
    if token[1] eq 'image' then infile = token[3]
    
    readf, iunit, line
    token = STRSPLIT( line, ' ', COUNT=count, /EXTRACT)
    if token[1] eq 'chip' then chip = fix(token[3])
    
    free_lun, iunit
    image = deimos_read_chip( infile, chip, header=header)
    
    ;; wait if requested...
    if i gt 0 and keyword_set(WAIT) then atv_activate

    ;; display...
    print, 'Displaying ', tracefile
    atv, image
    
    ;; read trace data...
    readcol, tracefile, u, v, comment='#', format='(f,f)'
    
    ;; show trace...
    atvplot, u, v, color='green'
    
endfor

if keyword_set(WAIT) then atv_activate

end

