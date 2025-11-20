;-----------------------------------------------------------------------
function strsub, template, search, replace
;-----------------------------------------------------------------------
; Purpose:
;       Replace the search string with the replace string and return
;
; Example:
;       psfile = 'foo.ps'
;       pdffile = strsub( psfile, 'ps', 'pdf')
;-----------------------------------------------------------------------
return, StrJoin( StrSplit( template, search, /Regex, /Extract, $
           /Preserve_Null), replace)
end
