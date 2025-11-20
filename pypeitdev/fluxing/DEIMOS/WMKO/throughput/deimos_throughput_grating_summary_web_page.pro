;-----------------------------------------------------------------------
pro deimos_throughput_grating_summary_web_page, outfile, grating, $
  params, summary
;-----------------------------------------------------------------------
;+
; NAME:
;	deimos_throughput_grating_summary_web_page
;
; PURPOSE:
;	
;       For each grating/grism, we generate the following page:
;       
;           * Summary plots
;                 o Current efficiency obtained by combining most
;                       recent dataset
;                 o Efficiency evolution show efficiency as a function
;                       of wavelength for various epochs
;                 o Efficiency evolution shows efficiency at a given 
;                       wavelength over time for various wavelengths 
;           * Data table
;                 o Setup
;                       + Grating
;                       + Central wavelength
;                       + Fitler 
;                 o Date
;                 o Star
;                 o Airmass
;                 o Measured efficiency at each passband
;                 o Link to detail page 
;       
;       Note: Table is sorted by date. 
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2009-Nov-15	GDW	Original version
;-

;; compile list of wavelengths for throughput checks...
cenlam = nint(params[0].lambda_eff)
n_cenlam = n_elements(cenlam)

;; define css...
css = 'TD.right { text-align: right; vertical-align: middle;}'

title = 'Throughput Data for Grating '+grating
body = ''

;;----------------------------------------
;; Figures...
;;----------------------------------------
body += '<h2>Current Efficiency</h2>'
body += 'Most recent plot of efficiency vs. wavelength (combine multiple spectra)'
body += '<p>'
body += '<a href="'+summary.efficiency_current_pdf+'">'
body += '<img src="'+summary.efficiency_current_plot+'" alt="efficiency plot">'
body += '</a>'

;;body += '<h2>Efficiency Evolution</h2>'
;;body += 'Efficiency as a function of wavelength for various epochs'

body += '<h2>Efficiency Evolution</h2>'
body += 'Efficiency as a function of time for various wavelengths'
body += '<p>'
body += '<a href="'+summary.eff_vs_time_pdf+'">'
body += '<img src="'+summary.eff_vs_time_plot+'" alt="efficiency plot">'
body += '</a>'

body += '<h2>Dataset Summary</h2>'
body += '<table border="1" cellpadding="5">'
body += '<thead>'
body += '<tr>'
body += '<th rowspan="2">Dataset</th>'
body += '<th rowspan="2">Date</th>'
body += '<th rowspan="2">Star</th>'
body += '<th rowspan="2">Airmass</th>'
body += '<th rowspan="2">Grating</th>'
body += '<th rowspan="2">Filter</th>'
body += '<th rowspan="2">&lambda;<sub>c</sub><br>[&Aring;]</th>'
body += '<th colspan="'+stringify(n_cenlam)+'">Efficiency(&lambda;) [%]</th>'
body += '</tr>'

;; label the efficiencies...
body += '<tr>'
    for j=0,n_cenlam-1 do begin
        body += '<th>'+stringify(cenlam[j])+'&nbsp;&Aring;</th>'
    endfor 
body += '</tr>'
body += '</thead>'
body += '<tbody>'

n = n_elements(params)
for i=0,n-1 do begin
    body += '<tr>'
    body += '<td><a href="'+params[i].detail+'">'+params[i].dataset+'</a></td>'
    body += '<td>'+params[i].date+'</td>'
    body += '<td>'+params[i].std_name+'</td>'
    body += '<td class="right">'+stringify(params[i].airmass,'(f15.2)')+'</td>'
    body += '<td>'+params[i].grating+'</td>'
    body += '<td>'+params[i].blocking+'</td>'
    body += '<td class="right">'+params[i].cenlam+'</td>'
    for j=0,n_cenlam-1 do begin
        if finite(params[i].efficiency[j]) then begin
            label = stringify(100.*params[i].efficiency[j],'(f15.2)')
        endif else begin
            label = '&nbsp;'
        endelse 
        body += '<td class="right">'+label+'</td>'
    endfor 
    body += '</tr>'
endfor 

body += '</tbody>'
body += '</table>'

;; write the file...
deimos_write_web_page, outfile, body, title=title, CSS=css

end
