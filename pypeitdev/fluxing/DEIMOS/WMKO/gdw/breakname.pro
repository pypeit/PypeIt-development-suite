;-----------------------------------------------------------------------
pro breakname, filename, dirname, rootname, extn
;-----------------------------------------------------------------------
;+
; NAME:
;	BREAKNAME
;
; PURPOSE:
; 	This procedure separates the the directory and extension names
; 	from a filename,  returning (1) the directory name, (2) the
; 	root name, (3) the extension.  The directory name ends just
; 	before the first character of the rootname; that is, at the
; 	last "/" character.  The extension begins at the first "."
; 	character after the rootname. 
;
; 	If a component is null, the empty string "" is returned.
;
; CATEGORY:
;	Filenames
;
; CALLING SEQUENCE:
;	BREAKNAME, Filename, Dirname, Rootname, Extn
;
; INPUTS:
;	Filename	[string] Name of a disk file to parse
;
; OUTPUTS:
;	Dirname		[string] Directory portion of the filename
;	Rootname	[string] Root portion of the filename
;	Extn		[string] Extension to the filename
;
; RESTRICTIONS:
;	Only works with Unix filenames
;
; EXAMPLE:
;	Extract the directory name, rootname, and extension from a
;	filename:
;
;		filename = '/usr/local/directory/root.extn'
;		breakname, filename, dirname, rootname, extn
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2000-Jan-20	GDW	Original version
;       2010-Sep-09     GDW     Rewrote to make everything after and
;       including the FIRST '.' in the filename be in the extension
;-
;-----------------------------------------------------------------------

; initialize...
l = strlen( filename)

; locate last character in directory string = i
i = l
found = 0B
for i=l-1,0,-1 do begin
    if strmid( filename, i, 1) eq '/' then begin
        found = 1B
        break
    endif 
endfor 
if found then begin
    dirname = strmid( filename, 0, i+1)
endif else begin
    dirname = ''
endelse 

; locate first character of rootname = j
found = 0B
for j=i+1,l-1 do begin
    if strmid( filename, j, 1) eq '.' then begin
        found = 1B
        break
    endif 
endfor 
if found then begin
    extn = strmid( filename, j)
endif else begin
    extn = ''
    j = l
endelse 

; extract rootname...
if i+1 lt j then begin
    rootname = strmid( filename, i+1, j-i-1)
endif else begin
    rootname = ''
endelse

end

