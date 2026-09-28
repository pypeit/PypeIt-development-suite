;-----------------------------------------------------------------------
pro deimos_throughput_grating_detail_web_page, params
;-----------------------------------------------------------------------
;+
; NAME:
;	DEIMOS_THROUHGPUT_GRATING_DETAIL_WEB_PAGE
;
; PURPOSE:
;	This procedure will generate a DEIMOS throughput grating
;	detail web page, providing the following information for a
;	single measurement:
;
;    *  Summary table (links for download)
;          o Efficiency vs. wavelength
;                + PNG
;                + PDF
;                + FITS bintab
;                + ASCII table 
;          o ZP(A) vs. wavelength
;                + PNG
;                + PDF
;                + FITS bintab
;                + ASCII table 
;          o ZP(pix) vs. wavelength
;                + PNG
;                + PDF
;                + FITS bintab
;                + ASCII table 
;    * Plots
;          o Efficiency vs. wavelength (with link to PNG of same)
;          o ZP (A) vs, wavelength (with link to PNG of same)
;          o ZP (px) vs, wavelength (with link to PNG of same) 
;
; CATEGORY:
;	Instrument performance monitoring
;
; CALLING SEQUENCE:
;       deimos_throughput_grating_detail_web_page, outfile, params
;
; INPUTS:
;	outfile:        name of the HTML file to generate (string)
;       params:         list of parameters (structure)
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
;       2009-Nov-15     GDW     Original version 
;-

title = 'Throughput Detail for Dataset '+params.dataset
body = ''
css = 'TD.right { text-align: right; vertical-align: middle;}'

;;----------------------------------------
;; list observing params...
;;----------------------------------------
body += '<h2>Observational Summary</h2>'
body += '<center>'
body += '<table cellpadding="5" border="1">'
body += '<thead>'
body += '<tr>'
body += '<th>Parameter</th>'
body += '<th>Value</th>'
body += '</tr>'
body += '</thead>'
body += '<tbody>'

body += '<tr>'
body += '<td>Dataset</td>'
body += '<td>'+params.dataset+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Date</td>'
body += '<td>'+params.date+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Star</td>'
body += '<td>'+params.std_name+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>RA</td>'
body += '<td>'+params.ra+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Dec</td>'
body += '<td>'+params.dec+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Grating</td>'
body += '<td>'+params.grating+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>&lambda;<sub>c</sub> [&Aring;]</td>'
body += '<td>'+params.cenlam+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Filter</td>'
body += '<td>'+params.blocking+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Slit Width [arcsec]</td>'
body += '<td>'+stringify(params.slit_width)+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Airmass</td>'
body += '<td>'+stringify(params.airmass)+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>SpecBin</td>'
body += '<td>'+params.spec_bin+'</td>'
body += '</tr>'

body += '<tr>'
body += '<td>Conditions</td>'
body += '<td>'+params.conditions+'</td>'
body += '</tr>'

body += '</tbody>'
body += '</table>'
body += '</center>'

;;----------------------------------------
;; list data products...
;;----------------------------------------
body += '<h2>Data Products</h2>'
body += '<center>'
body += '<table cellpadding="5" border="1">'
body += '<thead>'
body += '<tr>'
body += '<th>Format</th>'
body += '<th>Efficiency</th>'
body += '<th>Zero Point (&Aring;)</th>'
body += '<th>Zero Point (pix)</th>'
body += '</tr>'
body += '</thead>'
body += '<tbody>'

;; png...
body += '<tr>'
body += '<td>PNG Plot</td>'
body += '<td><a href="'+params.fig_eff+'">download PNG</a></td>'
body += '<td><a href="'+params.fig_zp_ang+'">download PNG</a></td>'
body += '<td><a href="'+params.fig_zp_pix+'">download PNG</a></td>'
body += '</tr>'

;; pdf...
body += '<tr>'
body += '<td>PDF Plot</td>'
body += '<td><a href="'+params.fig_eff_pdf+'">download PDF</a></td>'
body += '<td><a href="'+params.fig_zp_ang_pdf+'">download PDF</a></td>'
body += '<td><a href="'+params.fig_zp_pix_pdf+'">download PDF</a></td>'
body += '</tr>'

;; ascii...
body += '<tr>'
body += '<td>ASCII Table</td>'
body += '<td><a href="'+params.tab_eff+'">download ASCII</a></td>'
body += '<td><a href="'+params.tab_zp_ang+'">download ASCII</a></td>'
body += '<td><a href="'+params.tab_zp_pix+'">download ASCII</a></td>'
body += '</tr>'

;; ascii...
body += '<tr>'
body += '<td>FITS Binary Table</td>'
body += '<td colspan="3"><a href="'+params.infile+'">download FITS</a></td>'
body += '</tr>'

body += '</tbody>'
body += '</table>'
body += '</center>'

;;----------------------------------------
;; Figures...
;;----------------------------------------

;; compile list of wavelengths for throughput checks...
cenlam = params.lambda_eff
dcenlam = params.dlambda_eff
n_cenlam = n_elements(cenlam)
cenlam_label_min = stringify(nint(cenlam-0.5*dcenlam))
cenlam_label_max = stringify(nint(cenlam+0.5*dcenlam))
cenlam_label     = stringify(nint(cenlam))

body += '<h2>Efficiency</h2>'

body += '<a href="'+params.fig_eff_pdf+'">'
body += '<img src="'+params.fig_eff+'" alt="efficiency plot">'
body += '</a>'
body += '<p>'
body += '<center>'
body += '<table border="1" cellpadding="5">'
body += '<thead>'
body += '<tr>'
body += '<th>&lambda;<sub>c</sub><br>[&Aring;]</th>'
body += '<th>&lambda;<sub>min</sub><br>[&Aring;]</th>'
body += '<th>&lambda;<sub>max</sub><br>[&Aring;]</th>'
body += '<th>Efficiency<br>[%]</th>'
body += '</tr>'
body += '</thead>'
body += '<tbody>'

for j=0,n_cenlam-1 do begin
    body += '<tr>'
    body += '<td>'+cenlam_label[j]+'</td>'
    body += '<td>'+cenlam_label_min[j]+'</td>'
    body += '<td>'+cenlam_label_max[j]+'</td>'
    body += '<td class="right">'+stringify(100.*params.efficiency[j],'(f15.2)')+'</td>'
    body += '</tr>'
endfor 

body += '</tbody>'
body += '</table>'
body += '</center>'

body += '<h2>Zero Point [&Aring;]</h2>'
body += '<a href="'+params.fig_zp_ang_pdf+'">'
body += '<img src="'+params.fig_zp_ang+'" alt="efficiency plot">'
body += '</a>'

body += '<h2>Zero Point [pix]</h2>'
body += '<a href="'+params.fig_zp_pix_pdf+'">'
body += '<img src="'+params.fig_zp_pix+'" alt="efficiency plot">'
body += '</a>'

;; write the file...
deimos_write_web_page, params.detail, body, title=title, CSS=css

end
