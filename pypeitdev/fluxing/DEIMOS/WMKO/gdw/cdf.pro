;-----------------------------------------------------------------------
pro cdf, array, x, y
;-----------------------------------------------------------------------
; Given array of data values, generate and return a cumulative
; distribution function (CDF) which shows the fraction of array with
; values at or below the given X value.
;-----------------------------------------------------------------------

order = sort(array)
buf = array[order]

n = n_elements(array)
n2 = float(n)

x = fltarr(2*n)
y = fltarr(2*n)

j = 0
for i=0,n-1 do begin

    x[j] = buf[i]
    y[j] = float(i)/n2
    j += 1

    x[j] = buf[i]
    y[j] = float(i+1)/n2
    j += 1

endfor
end 
