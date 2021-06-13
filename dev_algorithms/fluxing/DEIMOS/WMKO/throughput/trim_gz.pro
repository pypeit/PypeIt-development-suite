;-----------------------------------------------------------------------
function trim_gz, filenames
;-----------------------------------------------------------------------
; purpose: trim .gz from filenames
;-----------------------------------------------------------------------

new = filenames
n = n_elements( new)
for i=0,n-1 do begin
    breakname, new[i], dirname, rootname, extn
    if extn eq '.fits.gz' then begin
        extn = '.fits'
        new[i] = dirname + rootname + extn
    endif 
    
endfor 

return, new
end
