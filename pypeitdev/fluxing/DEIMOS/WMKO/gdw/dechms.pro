;-----------------------------------------------------------------------
pro DECHMS, DH, PM, HH, MM, SS, LABEL=label, SIGN=sign
;-----------------------------------------------------------------------
;     Gregory D. Wirth / UC Santa Cruz / 19 Jan 1992
;     
;     Converts decimal hours to hours, minutes, seconds.
;
;      Implicit None
;      Real HH                   ! Hours
;      Real MM                   ! Minutes
;      Real SS                   ! Seconds
;      Real DH                   ! Decimal hours
;      Real D                    ! Absolute value of DH
;      Character*1 PM            ! Plus or minus character
;-----------------------------------------------------------------------

If DH lt 0. Then begin
    PM = '-'
endif Else begin
    PM = '+'
Endelse

D = Abs(DH)
HH = Double(Fix(D))
MM = Double(Fix((D-HH)*60.d0))
SS = (D - HH - MM/60.d0)*3600.d0

if dh lt 0. or keyword_set(SIGN) then begin
    label = string(format='(a,i02,":",i02,":",f04.1)', pm, hh, mm, ss)
endif else begin
    label = string(format='(i02,":",i02,":",f05.2)', hh, mm, ss)
endelse



End

