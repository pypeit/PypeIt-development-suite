;------------------------------------------------------------------------
function gdwtimestamp
;------------------------------------------------------------------------
; return a timestamp of the form YYYY-MM-DDTHH:MM:SS
;------------------------------------------------------------------------

buf = strsplit( systime(0), ' ', /extract)
year = buf[4]
mon = buf[1]
day = buf[2]
timestr = buf[3]

;; convert month to number...
case mon of
    'Jan': mm = 1
    'Feb': mm = 2
    'Mar': mm = 3
    'Apr': mm = 4
    'May': mm = 5
    'Jun': mm = 6
    'Jul': mm = 7
    'Aug': mm = 8
    'Sep': mm = 9
    'Oct': mm = 10
    'Nov': mm = 11
    'Dec': mm = 12
    else: message, mon+' is not a valid month'
ENDCASE
date = string( format='(i4,"-",i2.2,"-",i2.2,"T",a)', $
              year, mm, day, timestr )

return, date
end
