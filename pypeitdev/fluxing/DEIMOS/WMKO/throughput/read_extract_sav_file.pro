;------------------------------------------------------------------------
function read_extract_sav_file, file
;------------------------------------------------------------------------
; Purpose: given the name of a SAV file, read it and return the
; contents of the 'extract' struct.
;------------------------------------------------------------------------

;; create an initial value...
extract = 0B

;; the 'restore' command loads the contents of the SAV file.  If it
;; was correctly generated, it contains a single struct called 'extract'.
restore, file

if size( extract, /tname) ne 'STRUCT' then begin
    message, 'WARNING: there is no struct named EXTRACT in file '+file, /info
endif 

return, extract
end

