pro compare_throughput_files, in

n_in = n_elements(in)

;; first pass to determine wavelengths...
for i=0,n_in-1 do begin
    
    infile = in[i]
    spec = mrdfits( infile, 'SPECTRUM', header)
    
    ;; get range in wavelength and counts...
    x1 = min( spec.wavelength, max=x2)
    y1 = min( spec.counts, max=y2)

    if i eq 0 then begin
        w1 = x1
        w2 = x2
        c1 = y1
        c2 = y2
    endif else begin
        w1 = w1 < x1
        w2 = x2 > x2
        c1 = c1 < y1
        c2 = c2 > y2
    endelse 

endfor 

;; start plot...
psland, 'compare_throughput_files.ps'
DEVICE, /COLOR
TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
plot, [0], [0], xrange=[w1,w2], yrange=[c1,c2], /nodata

;; second pass to make plots...
for i=0,n_in-1 do begin
    
    infile = in[i]
    spec = mrdfits( infile, 'SPECTRUM', header)
    
    oplot, spec.wavelength, spec.counts, color=i
endfor 

;; clean up...
device, /close
set_plot, 'X'
end 
