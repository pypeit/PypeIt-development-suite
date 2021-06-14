;+
; NAME:
;	PSCONFIG_GDW
;
; PURPOSE:
;       This function provides a GUI interface to allow the user to
;       specify the destination for PostScript output, eerily similar
;       to that provided by the Netscape Navigator's "Print" widget.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	Result = PSCONFIG_GDW( GROUP_LEADER=group_leader )
;
; INPUTS:
;	None
;
; OPTIONAL INPUTS:
;       None
;	
; KEYWORD PARAMETERS:
;
;	GROUP_LEADER= (int)
;		Set this keyword to the ID number of a widget to serve
;		as the group leader.  For example, the top-level base
;		widget for the calling widget would be a suitable
;		choice. 
;
;	[ /TOPRINTER | /TOFILE ]
;		Set one of these keywords to select the default
;		destination.  The default is TOPRINTER.
;
;	[ /LANDSCAPE | /PORTRAIT ]
;		Set one of these keywords to select the default
;		orientation.  The default is LANDSCAPE.
;
;	/ENCAPSULATED
;		Set this keyword to change the default form of the
;		saved PostScript file to 'encapsulated'; the default
;		is to write non-encapsulated PostScript
;
;	[ /GREYSCALE | /COLOR ]
;		Set one of these keywords to set the default output
;		mode to either monochrome or color.  The default is
;		GREYSCALE. 
;
;	FILENAME= (string)
;		Set this keyword to the default name of the output
;		file.  The default is 'idl.ps'.
;
;	COMMAND= (string)
;		Set this keyword to the default system command to
;		use for processing output.  The default is 'enscript'.
;
;	XMARGIN= (real value; units of INCHES; default 1.0)
;	YMARGIN= (real value; units of INCHES; default 1.0)
;		Set either of these keywords to specify the page
;		margin, in inches.
;
;	[ /LETTER | /LEGAL | /EXECUTIVE | /A4 | /11x17]
;		Set one of these keywords to reset the default paper
;		size.  The default is LETTER.
;
; OUTPUTS:
;	This function returns a structure which contains the user's
;	requested PostScript output settings.  This structure can then 
;	be passed to the DEVICE procedure using the _EXTRA mechanism;
;	see example below.
;
; OPTIONAL OUTPUTS:
;
;	CANCELLED= (variable)
;		Set this keyword to the name of a variable which will
;		contain the value 1 if the user pressed the CANCEL
;		button, and 0 if the user exiting the GUI normally
;		(via the PRINT button).
;
; EXAMPLE:
;	The following code fragment demonstrates the use of this
;	procedure:
;		pskeywords = psconfig_gdw( cancelled=cancelled)
;		if cancelled then return
;	        if pskeywords.toprinter then begin
;	            buf = getenv("TMPDIR")
;	            if strlen(buf) eq 0 then buf = '/tmp'
;		    buf = buf + '/'
;	            pskeywords.filename = buf + 'idltemp' $
;	              + strtrim(string(long(systime(1))),2) $
;	              + ".ps"
;		thisDevice = !D.Name
;		set_plot, 'PS'
;		device, _EXTRA=pskeywords
;		 ...
;		plot, findgen(11) ; or other plot-generating commands...
;		 ...
;		device, /close
;		Set_Plot, thisDevice
;		if pskeywords.toprinter then begin
;		    command = pskeywords.command + ' ' + pskeywords.filename
;		    spawn, command, exit=status
;		    file_delete, pskeywords.filename
;		endif
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2002-Jan-11	GDW	Original version
;-

;-----------------------------------------------------------------------
pro psconfig_gdw_event, event
;-----------------------------------------------------------------------

widget_control, event.top, get_uvalue=ptr
widget_control, event.id, get_uvalue=uval
tlb=(*ptr).state.tlb

case uval of

    ; change mode...
    'mode': (*ptr).data.mode_index = event.value

    ; change mode...
    'paper': (*ptr).data.paper_index = event.value

    ; change output...
    'destination': begin

        ; store value...
        (*ptr).data.destination_index = event.value

        ; toggle sensitivity of other fields...
        printer_on = 1 - event.value
        file_on = event.value
        widget_control, (*ptr).state.command_entry, sensitive=printer_on
        widget_control, (*ptr).state.filename_entry, sensitive=file_on
        widget_control, (*ptr).state.filename_browse, sensitive=file_on
        widget_control, (*ptr).state.encapsulation_base, sensitive=file_on
        
    end

    'orientation': begin
        ; store value...
        (*ptr).data.orientation_index = event.value
    end

    'encapsulation': begin
        ; store value...
        (*ptr).data.encapsulation_index = event.value

        ; modify output filename...
        widget_control, (*ptr).state.filename_entry, get_value=buf
        filename = buf[0]
        if (*ptr).data.encapsulation_index eq 0 then begin
            old_suffix = ".eps"
            new_suffix = ".ps"
        endif else begin
            old_suffix = ".ps"
            new_suffix = ".eps"
        endelse

        i = strpos( filename, old_suffix, /REVERSE_SEARCH)
        j = strlen( old_suffix)
        k = strlen( filename)

        if i+j eq k then begin
            filename = strmid( filename, 0, i) + new_suffix
            widget_control, (*ptr).state.filename_entry, set_value=filename
        endif

    end

    ; browse...
    'browse': begin
        buf = dialog_pickfile( file=(*ptr).data.filename, $
                               /write, $
                               filter='*.ps')
        (*ptr).data.filename = buf
        widget_control, (*ptr).state.filename_entry, set_value=buf
    end

    ; print...
    'print': begin
        
        ; capture filename...
        widget_control, (*ptr).state.filename_entry, get_value=buf
        (*ptr).data.filename = buf

        ; capture command...
        widget_control, (*ptr).state.command_entry, get_value=buf
        (*ptr).data.command = buf

        widget_control, event.top, /destroy
    end

    ; cancel...
    'cancel': begin
        (*ptr).data.cancelled = 1
        widget_control, event.top, /destroy
    end

    ; no match...
    else:  message, uval + ' ' + string(event.id), /info

endcase

end

;-----------------------------------------------------------------------
function psconfig_gdw, TOPRINTER=toprinter, $
                       TOFILE=tofile, $
                       FILENAME=filename, $
                       COMMAND=command, $
                       GREYSCALE=greyscale, $
                       COLOR=color, $
                       LANDSCAPE=landscape, $
                       PORTRAIT=portrait, $
                       XMARGIN=xmargin, $
                       YMARGIN=ymargin, $
                       LETTER=letter, $
                       LEGAL=legal, $
                       EXECUTIVE=executive, $
                       A4=a4, $
                       LEDGER=ledger, $
                       CANCELLED=cancelled, $
                       GROUP_LEADER=group_leader, $
                       ENCAPSULATED=encapsulated
;-----------------------------------------------------------------------

; define constants...
label_width = 100

; define top level base...
tlb = widget_base( title='Print Menu', /col, /modal, $
                   group_leader=group_leader)

; define top section...
top = widget_base( tlb, /col, frame=1)

; define print destination...
destination_options = ['Printer','File']
destination_index = 0
if keyword_set(TOPRINTER) then destination_index = 0
if keyword_set(TOFILE)    then destination_index = 1
destination_base = widget_base( top, /row)
destination_label = WIDGET_LABEL( destination_base, $
                                  xsize=label_width, $
                                  VALUE='Print To:')
destination_bgroup = cw_bgroup( destination_base, $
                                destination_options, $
                                /no_release, $
                                /return_index, $
                                /exclusive, $
                                set_value=destination_index, $
                                uvalue='destination', $
                                /row)

; define print command...
sensitive = 1 - destination_index
if not keyword_set(COMMAND) then command = 'enscript'
command_base  = widget_base( top, /row)
command_label = WIDGET_LABEL( command_base, $
                              xsize=label_width, $
                              VALUE='Print Command:')
command_entry = widget_text( command_base, $
                             value=command, $
                             xsize=40, $
                             ysize=1, $
                             sensitive=sensitive, $
                             /editable)

; define file name...
sensitive = destination_index
if not keyword_set(FILENAME) then filename = 'idl.ps'
filename_base  = widget_base( top, /row)
filename_label = WIDGET_LABEL( filename_base, $
                               xsize=label_width, $
                               VALUE='File Name:')
filename_entry = widget_text( filename_base, $
                              value=filename, $
                              xsize=30, $
                              ysize=1, $
                              sensitive=sensitive, $
                              /editable)
filename_browse = widget_button( filename_base, $
                                 sensitive=sensitive, $
                                 value='Browse...', $
                                 uvalue='browse')

; define encapsulation mode...
sensitive = destination_index
encapsulation_options = ['Off','On']
encapsulation_index = 0
if keyword_set(ENCAPSULATED) then encapsulation_index = 1
encapsulation_base = widget_base( top, $
                                  /row, $
                                  sensitive=sensitive)
encapsulation_label = WIDGET_LABEL( encapsulation_base, $
                                    xsize=label_width, $
                                    VALUE='Encapsulation:')
encapsulation_bgroup = cw_bgroup( encapsulation_base, $
                                  encapsulation_options, $
                                  /no_release, $
                                  /return_index, $
                                  /exclusive, $
                                  set_value=encapsulation_index, $
                                  uvalue='encapsulation', $
                                  /row)

; define middle section...
middle = widget_base( tlb, /col, frame=1)

; define orientation button...
orientation_options = ['Portrait','Landscape']
orientation_index = 0
if keyword_set(PORTRAIT)  then orientation_index = 0
if keyword_set(LANDSCAPE) then orientation_index = 1
orientation_base = widget_base( middle, /row)
orientation_label = WIDGET_LABEL( orientation_base, $
                                  xsize=label_width, $
                                  VALUE='Orientation:')
orientation_bgroup = cw_bgroup( orientation_base, $
                                orientation_options, $
                                /no_release, $
                                /return_index, $
                                /exclusive, $
                                set_value=orientation_index, $
                                uvalue='orientation', $
                                /row)

; define mode button...
mode_options = ['Greyscale','Color']
mode_index = 0
if keyword_set(GREYSCALE) then mode_index = 0
if keyword_set(COLOR)     then mode_index = 1
mode_base = widget_base( middle, /row)
mode_label = WIDGET_LABEL( mode_base, $
                           xsize=label_width, $
                           VALUE='Mode:')
mode_bgroup = cw_bgroup( mode_base, $
                         mode_options, $
                         /no_release, $
                         /return_index, $
                         /exclusive, $
                         set_value=mode_index, $
                         uvalue='mode', $
                         /row)

; define paper button...
paper_options = ['Letter (8 1/2 x 11 in.)', $
                 'Legal (8 1/2 x 14 in.)', $
                 'Executive (7 1/2 x 10 in.)', $
                 'Ledger (11 x 17 in.)', $
                 'A4 (210 x 297)' ]
shortside_options = [8.5, 8.5, 7.5, 11, 8.27]
longside_options = [11, 14, 10, 17, 11.7]
paper_index = 0
if keyword_set(LETTER)    then paper_index = 0
if keyword_set(LEGAL)     then paper_index = 1
if keyword_set(EXECUTIVE) then paper_index = 2
if keyword_set(LEDGER)    then paper_index = 3
if keyword_set(A4)        then paper_index = 4
paper_base = widget_base( middle, /row)
paper_label = WIDGET_LABEL( paper_base, $
                           xsize=label_width, $
                           VALUE='Paper Size:')
paper_bgroup = cw_bgroup( paper_base, $
                          paper_options, $
                          /no_release, $
                          /return_index, $
                          /exclusive, $
                          set_value=paper_index, $
                          uvalue='paper')

; define print and cancel buttons...
button_base = widget_base( tlb, /row)
print_button = widget_button( button_base, value='Print', $
                                 uvalue='print')
cancel_button = widget_button( button_base, value='Cancel', $
                                 uvalue='cancel')

; make the print button the default...
widget_control, tlb, $
  default_button=print_button, $
  cancel_button=cancel_button

; define state structure...
state = { tlb:tlb, $
          command_entry:command_entry, $
          filename_entry:filename_entry, $
          filename_browse:filename_browse, $
          encapsulation_base:encapsulation_base}

data = { destination_index:destination_index, $
         orientation_index:orientation_index, $
         encapsulation_index:encapsulation_index, $
         mode_index:mode_index, $
         paper_index:paper_index, $
         filename:filename, $
         command:command, $
         cancelled: 0}

; define global structure...
info = { state:state, $
         data:data }

; Create a pointer to the info structure.  Instead of passing around
; the entire structure, this is all that will be passed to the event
; handlers...
ptr = ptr_new( info, /no_copy)
widget_control, tlb, /realize, set_uvalue=ptr

; Register the widget with Xmanager and pause execution until
; user completes input on the GUI...
xmanager, 'psconfig_gdw', tlb

; unpack data...
toprinter = 1 - (*ptr).data.destination_index
tofile = (*ptr).data.destination_index

landscape = (*ptr).data.orientation_index
portrait = 1 - (*ptr).data.orientation_index

greyscale = 1 - (*ptr).data.mode_index
color = (*ptr).data.mode_index

shortside = shortside_options((*ptr).data.paper_index)
longside = longside_options((*ptr).data.paper_index)

if tofile then begin
    encapsulated = (*ptr).data.encapsulation_index
endif else begin
    encapsulated = 0
endelse

; define margins...
if not keyword_set(XMARGIN) then xmargin = 1.
if not keyword_set(YMARGIN) then ymargin = 1.

; flip page sizes if landscape...
if portrait then begin
    xsize = shortside
    ysize = longside
    xoffset = xmargin
    yoffset = ymargin
endif else begin
    xsize = longside
    ysize = shortside
    xoffset = ymargin
    yoffset = longside - xmargin
end

; reduce plot size by margins...
xsize = xsize - 2*xmargin
ysize = ysize - 2*ymargin

; if the CANCEL button was pressed, set the appropriate output flag...
cancelled = (*ptr).data.cancelled

; create struct...
struct = { toprinter: toprinter, $
           tofile: tofile, $
           filename: (*ptr).data.filename, $
           command: (*ptr).data.command, $
           encapsulated: encapsulated, $
           greyscale: greyscale, $
           color: color, $
           landscape: landscape, $
           portrait: portrait, $
           xsize: xsize, $
           ysize: ysize, $
           xoffset: xoffset, $
           yoffset: yoffset, $
           inches: 1 }

return, struct

end
