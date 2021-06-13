pro do_compare_throughput_files

dir = '/s/sdata1001/dmoseng/kroot/throughput/extract/830G/'
in = ['2008oct02_d1002_0049.fits','2008oct02_d1002_0050.fits']
in = dir + in
compare_throughput_files, in

;dir = '/h/scratch13/gwirth/deimos-throughput/'
;in = ['spec0049.fits','spec0050.fits']
;in = dir + in
;compare_throughput_files, in

end

