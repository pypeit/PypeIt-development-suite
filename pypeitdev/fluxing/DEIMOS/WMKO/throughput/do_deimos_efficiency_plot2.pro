thru_dir = getenv('DEIMOS_THRU_DIR')

indir  = thru_dir + '/eff/600ZD'
outdir = thru_dir + '/result/600ZD'
infiles = FILE_SEARCH(indir+'/*.sav')

graname = '600ZD'
;; cenlams = [5500, 6500, 7500]
lambda_eff = [5000., 6200., 8000., 9000.]

;; loop over lambda_eff values...
nl = n_elements(lambda_eff)
for i=0,nl-1 do begin
    
    l = nint(lambda_eff[i])
    format = '(a,"_",i4,".eps")'
    psfile = outdir + '/' $
             + string( format=format, graname, l)
    
    deimos_efficiency_plot, infiles, psfile, $
      GRANAME=graname, LAMBDA_EFF=l
endfor 

end
