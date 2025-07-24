function xDateString2JD, xString
;; purpose: convert a Prochaska-style date string to JD

;; perform pattern match...
ucString = strupcase(xString)
pattern = '^([0-9][0-9])([A-Z][A-Z][A-Z])([0-9][0-9][0-9][0-9])$'
buf = stregex( ucString, pattern, /subexpr, /extract )

;; check for no match...
if strlen(buf[0]) lt 1 then begin
    message, 'invalid date string', /info
    return, !values.f_nan
endif 
    
;; grab components...
day   = buf[1]
month = month_cnv(buf[2])       ; integer
year  = buf[3]

;; convert to string...
date = string( format='(i4,"-",i2.2,"-",i2.2)', year, month, day)

;; convert to JD...
jd = date_conv(date, 'J')

;; convert back to year...
;date2 = date_conv( jd, 'F')
;print, xString, ' ', date, ' ', jd, ' ', date2

return, jd

end
