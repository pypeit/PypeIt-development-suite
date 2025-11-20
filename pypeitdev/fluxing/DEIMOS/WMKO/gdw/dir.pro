pro dir,filename
;+
; NAME:
;	DIR
; PURPOSE:
;	Emulates the DCL command DIRECTORY on a VMS machine
;	(On a UNIX cmachine DIR will simply spawn an ls command.)
; CALLING SEQUENCE:
;	DIR  [, FILENAME ]
; OPTIONAL INPUT PARAMETERS:
;	FILENAME - Filename specification.  Can be of any form that DCL or
;                  UNIX recognizes, with some provisos listed below.
;                  If omitted, then "." is assumed for Unix, or
;                  "*.*;*" is assumed if on a VAX.  No parsing of 
;                  the filename is performed in Unix.
;                                         
;	The following conventions should be used:
;
;	Use DIR,'*.FOR' instead of DIR,'.FOR' to list all files with the 
;	extension ".FOR".
;       
;	Use DIR,'IDL$SYSPROG:' instead of DIR,'IDL$SYSPROG' to list all 
;	files in the directory IDL$SYSPROG on a VAX.
;
; OUTPUTS:
;	None.
; UNIX PROCEDURE:
;	The filename specification is formed, if necessary, and a command is
;	spawned to the operating system to list the directory.
; VAX PROCEDURE:
;	The filename specification is formed, and FINDFILE is used to find the
;	filenames.  These are then parsed to produce output at the terminal.
; MODIFICATION HISTORY:
;	William Thompson	Applied Research Corporation
;	September, 1987		8201 Corporate Drive
;				Landover, MD  20785
;	Changed for use on a SUN workstation from the VMS version.  
;		M. Greason, STX, August 1990.
;       Combined VAX and SUN versions, N. Collins, Nov 1990.
;-
;
if !VERSION.ARCH ne "vax" then begin
;
;  SUN version
;
;  			If no filename was passed, assume ".".
;
	IF N_PARAMS(0) EQ 1 THEN FILE = FILENAME ELSE FILE = "."
;
;			Spawn to the operating system.
;
	spawn, 'ls ' + file
;

endif else begin                    ; End SUN version
;
;  VAX version
;
;  If no filename was passed, assume "*.*;*".
;
	IF N_PARAMS(0) EQ 1 THEN FILE = FILENAME ELSE FILE = "*.*;*"
;
;  Strip off any preceeding directory, recognized by a bracket "]" or by a
;  colon ":", and look to see if the filename contains a file extension, 
;  recognized by a period.  If no file extension has been passed, assume an 
;  extension of ".*".  If only the directory name has been passed, list all
;  the files in that directory.
;
	START = (STRPOS(FILE,"]") > STRPOS(FILE,":")) + 1
	IF START EQ STRLEN(FILE) THEN FILE = FILE + "*.*;*" ELSE BEGIN
		PERIOD = STRPOS(FILE,".",START)
		IF PERIOD EQ -1 THEN FILE = FILE + ".*"
	ENDELSE
;
;  Look for any version numbers, recognized by a semi-colon ";".  If none was 
;  passed, assume ";*"
;
	IF STRPOS(FILE,";") EQ -1 THEN FILE = FILE + ";*"
;
;  Find all the filenames answering that description.
;
	FILES = FINDFILE(FILE,COUNT=N_FILES)
        IF N_FILES EQ 0 THEN BEGIN
              PRINT,'DIR - No files found - ',spec_dir(FILE)
              RETURN
        END
;
;  Set up all the variables necessary to begin parsing the filenames.
;
	NAMES = STRARR(4)               ;Array to store filenames in.
	J_FILES = 0			;Number of filenames in current line.
	N_SECTIONS = 0			;Number of A20 sections currently used.
	N_DIRECTORIES = 0		;Number of directories counted.
	TOTAL_FILES = 0			;Number of files in current directory.
	LAST_DIRECTORY = ''		;Previous directory name
	FORMAT = '(1X'			;Partial format specification.
;
;  Check each returned filename in order.  Extract from that filename the 
;  directory name and the rest of the filename.
;
	FOR I_FILE = 0,N_FILES-1 DO BEGIN
		FILE = FILES(I_FILE)
		BRACKET = STRPOS(FILE,"]") + 1
		DIRECTORY = STRMID(FILE,0,BRACKET)
		NAME = STRMID(FILE,BRACKET,STRLEN(FILE)-BRACKET)
;
;  If the directory name of the file is different from the previous directory 
;  name, print out any information left over from the last directory, and print 
;  the header for the new directory.
;
		IF DIRECTORY NE LAST_DIRECTORY THEN BEGIN
			IF TOTAL_FILES NE 0 THEN BEGIN
;
;  Print any remaining files.
;
				IF N_SECTIONS NE 0 THEN BEGIN
					FORMAT = FORMAT + ')'
		                        PRINT,FORMAT=FORMAT,NAMES(0:J_FILES-1)
					FORMAT = '(1X'
					J_FILES = 0
					N_SECTIONS = 0
				ENDIF
;
;  Increment the directory number, and print the trailer.
;
			    N_DIRECTORIES = N_DIRECTORIES + 1
			    PRINT,' '
			    PRINT,'Total of ',STRTRIM(TOTAL_FILES,2),' files.'
			ENDIF
;
;  Print the header.
;
			PRINT,' '
			PRINT,'Directory ',DIRECTORY
			PRINT,' '
			TOTAL_FILES = 0
			LAST_DIRECTORY = DIRECTORY
		ENDIF
;
;  Check to see how many characters the filename would need, in units of 20.
;
		I_SECTION = ((FIX(STRLEN(STRTRIM(NAME,2)))/20) + 1) < 4
;
;  If that would overflow the current line, print the line and start a new one.
;
		IF (N_SECTIONS + I_SECTION) GT 4 THEN BEGIN
			FORMAT = FORMAT + ')'
			PRINT,FORMAT=FORMAT,NAMES(0:J_FILES-1)
			FORMAT = '(1X'
			J_FILES = 0
			N_SECTIONS = 0
		ENDIF
;
;  Add the filename to the current line, and augment the format specification 
;  accordingly.  For example, if the filename contains less than 20 characters,
;  add the format specification ',A20' to the variable FORMAT.
;
		NAMES(J_FILES) = NAME
		J_FILES = J_FILES + 1
		N_SECTIONS = N_SECTIONS + I_SECTION
		TOTAL_FILES = TOTAL_FILES + 1
		FORMAT = FORMAT + ',A,T' + STRTRIM(20*N_SECTIONS+1,2)
	ENDFOR
;
;  Print out the remaining line, and the last trailer.
;
	FORMAT = FORMAT + ')'
	PRINT,FORMAT=FORMAT,NAMES(0:J_FILES-1)
	PRINT,' '
	PRINT,'Total of ',STRTRIM(TOTAL_FILES,2),' files.'
;
;  If more than one directory was listed, print out the grand total.
;
	N_DIRECTORIES = N_DIRECTORIES + 1
	IF N_DIRECTORIES GT 1 THEN BEGIN
		PRINT,' '
	PRINT,'Grand total of ',STRTRIM(N_DIRECTORIES,2),  $
	      ' directories, ',STRTRIM(N_FILES,2),' files.'
	ENDIF
endelse
;
	RETURN
	END
