;-----------------------------------------------------------------------
pro build_throughput_file_list_html
;-----------------------------------------------------------------------
;+
; NAME:
;	BUILD_THROUGHPUT_FILE_LIST_HTML
;
; PURPOSE:
;	Reads the input.lst and and generates input.html, a more
;	readable format sorted by grating and target.
;
; CATEGORY:
;	IPM
;
; CALLING SEQUENCE:
;	build_throughput_file_list_html
;
; INPUTS:
;       None
;
; OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLE:
;
; NOTES:
;       After writing this, I decided the better idea was to read
;       input.lst into excel.
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2012-Jan-07	GDW	Original version
;-
;-----------------------------------------------------------------------

;; define files...
thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq '' then $
  message, 'DEIMOS_THRU_DIR is not defined -- abort!'
data_dir = thru_dir + '/dat/'
input = data_dir + 'input.lst'
output = data_dir + 'output.html'

;; thou shalt not kill...
if file_test(output) then $
  message, 'ERROR: operation would clobber existing file '+output

;; read input file...
buf = read_csv( input)

;; re-organize struct...
n = n_elements( buf.field1)
s = {filename:'', targname:'', gratenam:'', wavelen:0L, dwfilnam:'', $
     dateobs:'', status:0L }
info = replicate(s, n)

for i=0,n-1 do begin
    info[i] = { filename:buf.field1[i], $
                targname:buf.field2[i], $
                gratenam:buf.field3[i], $
                wavelen: buf.field4[i], $
                dwfilnam:buf.field5[i], $
                dateobs: buf.field6[i], $
                status:  buf.field7[i] }
endfor 

;; sort fields by grating, target, and date...
order = sort( info.targname)
info = info[order]

order = sort( info.gratenam)
info = info[order]

;; write results...
openw, ounit, output, /get_lun

;; init file...
printf, ounit, '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">'
printf, ounit, '<html>'
printf, ounit, '  <head>'
printf, ounit, '    <title>DEIMOS Throughput Summary</title>'
printf, ounit, '  </head>'
printf, ounit, '  <body>'
printf, ounit, '    <h1>DEIMOS Throughput Summary</h1>'
printf, ounit, '    <table border="1">'
printf, ounit, '    <thead>'
printf, ounit, '    <tr>'
printf, ounit, '    <th>No.</th>'
printf, ounit, '    <th>Filename</th>'
printf, ounit, '    <th>Targname</th>'
printf, ounit, '    <th>Gratenam</th>'
printf, ounit, '    <th>Wavelen</th>'
printf, ounit, '    <th>Dwfilnam</th>'
printf, ounit, '    <th>DateObs</th>'
printf, ounit, '    <th>Status</th>'
printf, ounit, '    </tr>'
printf, ounit, '    </thead>'
printf, ounit, '    <tbody>'

;; loop over entries...
for i=0,n-1 do begin
    l = info[i]
    printf, ounit, '<tr>'
    printf, ounit, '<td>', i, '</td>'
    printf, ounit, '<td>', l.filename, '</td>'
    printf, ounit, '<td>', l.targname, '</td>'
    printf, ounit, '<td>', l.gratenam, '</td>'
    printf, ounit, '<td>', l.wavelen, '</td>'
    printf, ounit, '<td>', l.dwfilnam, '</td>'
    printf, ounit, '<td>', l.dateobs, '</td>'
    printf, ounit, '<td>', l.status, '</td>'
    printf, ounit, '</tr>'
endfor 

;; complete file...
printf, ounit, '    </tbody>'
printf, ounit, '    </table>'
printf, ounit, '    <p><hr>'
printf, ounit, 'Created on ', systime(0)
printf, ounit, '</body>'
printf, ounit, '</html>'

;; close output...
free_lun, ounit

end
