;-----------------------------------------------------------------------
function scan, buf, word
;-----------------------------------------------------------------------
; GDW / WMKO / 15 Sep 1998
;
; Purpose:
;     Return the next word from a character buffer, and remove the
;     word from the buffer.  Returned value is the number characters
;     in the word.
;
; Usage:
;     nchar = scan( string, word)
;-----------------------------------------------------------------------

; remove leading spaces...
while strpos( buf, ' ') eq 0 do buf = strmid( buf, 1)

; search for next space...
i = strpos( buf, ' ')
if( i lt 0) then begin
    word = buf
    buf = ''
endif else begin
    word = strmid( buf, 0, i)
    buf = strmid( buf, i+1)
endelse

; return the number of characters scanned...
return, strlen( word)

end
