;-----------------------------------------------------------------------
function raw2extract, infile, HEADER=header, SUBDIR=subdir
;-----------------------------------------------------------------------

;; determine output subdirectory...
outdir = getenv('DEIMOS_THRU_DIR')

;; read header as needed..
if ~ keyword_set(HEADER) then header = headfits( infile)

;; determine grating...
gratenam = strtrim(sxpar( header, 'gratenam'),2)

;; construct directory name...
subdir = outdir + '/extract/' + gratenam

;; build output filename based on input dirname...
dirname = file_basename(file_dirname(infile))
basename = file_basename(infile)
buf = strsplit( basename, '.', /extract)
root = buf[0]
stem = subdir + '/' + dirname + '_' + root

return, stem

end
