pro do_deimos_throughput_faber_plots

; cd, '/s/sdata1001/dmoseng/kroot/throughput/extract/plots/'
dir = '/s/sdata1001/dmoseng/kroot/throughput/extract/1200G/'

list = 'BD+28deg4211_1200G_slitless_pa0.extract.lst'
readcol, list, infiles, format='(a)', comment='#'
infiles = dir + infiles

;;list = 'BD+28deg4211_1200G_slitless_pa0.extract.psfiles.lst'
;;readcol, list, psfiles, format='(a)', comment='#'
psfiles = infiles
for i=0,n_elements(psfiles)-1 do begin
    breakname, infiles[i], dir, root, extn
    psfiles[i] = root + '.ps'
endfor 

deimos_throughput_faber_plots, infiles, psfiles, 5000
deimos_throughput_faber_plots, infiles, psfiles, 7000
deimos_throughput_faber_plots, infiles, psfiles, 8000
end
