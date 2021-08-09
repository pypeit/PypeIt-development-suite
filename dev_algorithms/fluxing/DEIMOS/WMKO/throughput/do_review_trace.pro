;; infile = '/s/sdata1001/dmoseng/kroot/throughput/extract/1200G/files'
infile = '/s/sdata1001/dmoseng/kroot/throughput/misc/BD+28deg4211_1200G_slitless_pa0.regions'
readcol, infile, filelist, format='(a)', comment='#'
review_trace, filelist, /wait
end
