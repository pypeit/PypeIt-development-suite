;-----------------------------------------------------------------------
function sex2dec, array
;-----------------------------------------------------------------------
; Converts from sexagesimal (e.g. +HH:MM:SS.SS) to decimal
;-----------------------------------------------------------------------

n = n_elements(array)
result = dblarr(n)
for i=0,n-1 do begin
    s = array[i]
    if (strmid(s,0,1) eq "-") then pm=-1. else pm=1.
    parts = strsplit( s, ':', /extract)
    m = n_elements(parts)
    sum = 0.d0
    for j=0,m-1 do begin
        sum = sum + abs(float(parts[j])/60.^j)
    endfor
    result[i] = pm*sum
endfor

return, result

end
