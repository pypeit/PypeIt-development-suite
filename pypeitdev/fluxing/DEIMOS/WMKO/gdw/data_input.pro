;-----------------------------------------------------------------------
pro data_input, datafile, maxlines, numcolumns, numlines, data, $
                VERBOSE=verbose, COMMENT=comment
;-----------------------------------------------------------------------
; Gregory D. Wirth / Northwestern University / July 7 1987
;
; This IDL procedure produces arrays of string data that correspond to the
; columns in the specified input file.  Any header lines beginning
; with the comment character ('#' by default) are skipped.
; If requested with /verbose, the procedure also prints
; out the total number of lines of data and the file from which
; they are read.
;
; Parameters:
; 	string	datafile	I: name of the input data file
; 	int	maxlines	I: maximum number of elements in each array
; 	int	NumColumns	I: number of columns of data in file Datafile
; 	int	NumLines	O: number of lines of data read from file
; 	string	Data[numcolumns,numlines] O: data array
;
; Usage:
;       data_input, filename, maxlines, numcolumns, numlines, data
;       a = reform(data[0,*])
;       b = float( reform(data[1,*]))
;       c = fix( reform(data[1,*]))
;       ...
;-----------------------------------------------------------------------
verbose = keyword_set( VERBOSE)

; define the comment character...
if( keyword_set( COMMENT)) then $
  comment_char = comment $
else $
  comment_char = '#'

; open file...
openr, unit, datafile, /get_lun

; initialize buffers as type string...
line = ''
data = strarr(numcolumns,maxlines)

; initialize line counter...
numlines = 0

; loop over lines in file...
while not eof(unit) and numlines le maxlines do begin
    readf, unit, line

    ; skip comment lines...
    if( strpos( line, comment_char) ne 0) then begin
        numlines = numlines + 1

        ; extract words from data...
        for i=0,numcolumns-1 do begin	   
            istat = scan( line, word)
            data(i, numlines-1) = word
        endfor

    endif
endwhile

; close file...
free_lun, unit

; truncate the unused elements...
data = data(*,0:numlines-1)

; print message...
if verbose then print, numlines,' lines read from ', datafile
 
end
