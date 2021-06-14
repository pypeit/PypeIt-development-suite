pro do_deimos_sensstd

indir = '/s/sdata1001/dmoseng/throughput/extract'
infiles = file_search( indir, '*fits*')
n_in = n_elements(infiles)

;; loop over input files...
for i=0,n_in-1 do begin
    in = infiles[i]
    message, 'processing '+in, /info
    deimos_sensstd_gdw, in, /no_upd_web
endfor 

print, 'processed ', n_in, ' files'

end
