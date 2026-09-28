;-----------------------------------------------------------------------
function dec2sex, value, SHORT=short, NOSIGN=nosign
;-----------------------------------------------------------------------
; Convert decimal value to sexigesimal.
;-----------------------------------------------------------------------

if value lt 0 then begin
    pm = "-"
endif else begin
    pm = "+"
endelse

a = abs(value)
dd = fix(a)
mm = fix((a - dd)*60.)
ss = (a - dd - mm/60.)*3600.
if keyword_set(SHORT) then begin
    result = string(format='(i02,":",i02)',dd,nint(mm+ss/60.))
endif else begin
    result = string(format='(i02,":",i02,":",f05.2)',dd,mm,ss)
endelse
if not keyword_set(NOSIGN) then result = pm + result
return, result
end
