;-----------------------------------------------------------------------
function mktemp, filename
;-----------------------------------------------------------------------
; Purpose:
;       When passed a template filename string with a '%' in it,
;       substitute a timestamped string and return the new filename.
;       Used to create a nearly unique temp file name.
;
; Example:
;       IDL> foo = 'test%.ps'
;       IDL> filename = mktemp(foo)
;       IDL> print, filename
;       test1256954587322.ps
;------------------------------------------------------------------------

substring = strtrim(string(ulong64(systime(/sec)*1000.)),2)
return, StrJoin( StrSplit( filename, '\%', /Regex, /Extract, $
           /Preserve_Null), substring)

end
