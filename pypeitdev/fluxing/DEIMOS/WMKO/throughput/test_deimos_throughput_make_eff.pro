thru_dir = getenv('DEIMOS_THRU_DIR')

indir  = thru_dir + '/calibs/600ZD'
outdir = thru_dir + '/eff/600ZD'

infiles = FILE_SEARCH(indir+'/*.sav')
n = n_elements( infiles)

file_mkdir, outdir

for i=0,n-1 do begin
    in = infiles[i]
    breakname, in, dir, root, extn
    out = outdir + '/' + root + extn

    ;; print, in, ' > ', out
    deimos_throughput_make_eff, in, out

endfor 

end

