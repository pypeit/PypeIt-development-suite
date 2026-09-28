pro make_summary_database, outfile
;+
; NAME:
;	make_summary_database
;
; PURPOSE:
;	Generate a FITS table summarizing the input spectra.
;
; CATEGORY:
;	IPM
;
; CALLING SEQUENCE:
;       make_summary_database, outfile
;
; INPUTS:
;	outfile: name of the output fits table to create
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2012-Jan-07	GDW	Original version
;-

;; get file list...
indir = '/s/sdata1001/dmoseng/kroot/throughput/raw'
files = FILE_SEARCH( indir+'/*/*.fits*')

;; generate a struct...
a = {filename:'', gratenam:'', wavelen:'', dwfilnam:'', targnam:'', dateobs:''}
data = replicate( a, n_files)

;; loop over files
n = n_elements(files)
for i=0,n-1 do begin

    ;; read header...
    f = files[i]
    h = mrdfits( f, 0, /silent)

    ;; print results...
    print, f, h.grating, nint(h.central_wave), h.blocking, h.std_name, format=format
    printf, ounit, f, h.grating, nint(h.central_wave), h.blocking, h.std_name, format=format
endfor 

end
