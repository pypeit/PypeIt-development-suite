;-----------------------------------------------------------------------
pro deimos_write_web_page, outfile, body, TITLE=title, INDEX=index, CSS=css
;-----------------------------------------------------------------------
;+
; NAME:
;	DEIMOS_WRITE_WEB_PAGE
;
; PURPOSE:
;       This procedure will write the passed body into a DEIMOS-format
;       web page, adding appropriate header and footer information.
;
; CATEGORY:
;	Web doco
;
; CALLING SEQUENCE:
;       DEIMOS_WRITE_WEB_PAGE, outfile, title, index, body
;
; INPUTS:
;	outfile:        Name of the HTML file to write.
;
;       title:          Title of the web page
;
;       index:          Index entry for the web page
;
;       body:           Body of the web page
;
; EXAMPLE:
;       1) Write the HTML code in the variable "body" to the file
;       throughput.html:
;               deimos_write_web_page, 'throughput.html', $
;                       'Throughput Measurements', $
;                       
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2000-Jan-01	GDW	Original version
;-

;-----------------------------------------------------------------------

;; build optional index string...
if keyword_set(INDEX) then begin
  index_string = '<!-- index="'+index+'" -->' 
endif else begin
  index_string = ''
endelse 

if keyword_set(TITLE) then begin
    title_string = title 
endif else begin
    title_string = ''
endelse 

;; Header
header = '<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">'
header += '<html><head><title>' + title_string + '</title></head>'
header += '<style>'
header += 'body  {background:cornsilk;}'
header += 'table {background-color:white;}'
header += 'thead {background-color:wheat;}'

;; add in CSS...
if keyword_set(CSS) then begin
    header += css
endif

header += '</style>'
header += ''
header += index_string
header += '<body>'
header += ''
header += '<table border="0" width="100%" >'
header += '<tr>'
header += '<td width="1%" bgcolor="#000000" align="center">'
header += '<img src="http://www.keck.hawaii.edu/realpublic/inst/deimos/logo.jpg" align="center">'
header += '</td>'
header += ''
header += '<td bgcolor="#000000">'
header += '<center>'
header += '<b><font color="#ffffff" size="+3">DEIMOS</font></b>'
header += '<br>'
header += '<font size="+2" color="#ffffff">' + title_string + '</font>'
header += '</center>'
header += '</td>'
header += ''
header += '<td width="1%" bgcolor="#000000" align="center">'
header += '<img src="http://www.keck.hawaii.edu/realpublic/inst/deimos/logo.jpg" align="center">'
header += '</td>'
header += '</tr>'
header += '</table>'
header += '<div style="margin-left:70px;margin-right:70px;">'

footer = '<p><hr>' 
footer += '<center>' 
footer += '<nobr>Go to:</nobr>'
footer += '<nobr><a href="index.html">DEIMOS Home Page</a></nobr> -' 
footer += '<nobr><a href="../">Instruments Home Page</a></nobr> -' 
footer += '<nobr><a href="../../..">Keck Home Page</a></nobr>' 
footer += '</center>' 
footer += '<hr>' 
footer += '<address>' 
footer += '<script language=javascript>' 
footer += '<!--' 
footer += 'var contact = "Keck Instrument Group"' 
footer += 'var email = "instrument"' 
footer += 'var emailHost = "keck.hawaii.edu"' 
footer += 'document.write("<a href=" + "mail" + "to:" + email + "@" + emailHost+ ">" + contact + "</a>")' 
footer += '//-->' 
footer += '</script>' 
footer += '</address>' 
footer += '<!-- hhmts start -->' 
footer += 'Last modified: '+systime(/utc)
footer += '<!-- hhmts end -->' 
footer += '</div>' 
footer += '</body>' 
footer += '</html>'

;; write to file...
openw, outunit, outfile, /get_lun
printf, outunit, header, body, footer
free_lun, outunit

end
