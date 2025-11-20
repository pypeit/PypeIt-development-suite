;-----------------------------------------------------------------------
function get_wavelength_shift, wavelength1, counts1, $
  wavelength2, counts2, status=status, chisq=chisq
;-----------------------------------------------------------------------
;+
; NAME:
;	get_wavelength_shift
;
; PURPOSE:
;	This function will compare the passed spectrum with the
;	fiducial and return the wavelength shift which brings them
;	into alignment.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;       test_wavelength_shift, wavelength1, counts1, wavelength2, counts2
;
; INPUTS:
;	wavelength1: array of wavelenths for the fiducial spectrum
;	counts1:     array of counts (signal) for the fiducial spectrum
;	wavelength1: array of wavelenths for the fiducial spectrum
;	counts1:     array of counts (signal) for the fiducial spectrum
;
; OUTPUT KEYWORDS:
;       status: status of the fit (0=good, non-zero=bad)
;       chisq: goodness of fit parameter
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; AUTHOR:
;	
;
; MODIFICATION HISTORY:
; 	
;-
;-----------------------------------------------------------------------

status = 0
delta = -100.
return, delta

end

;-----------------------------------------------------------------------
pro doit, file1, file2, SHIFT=shift
;-----------------------------------------------------------------------
;+
; NAME:
;	doit
;
; PURPOSE:
;	This function will read the spectra in the two specified 
;       filenames, plot them, invoke the get_wavelength_shift function
;       to compute the offset, and then will plot the result.
;
; CALLING SEQUENCE:
;      doit, file1, file2
;
; INPUTS:
;	file1:	Name of the FITS image containing the reference spectrum
;	file2:	Name of the FITS image containing the spectrum to shift
;
; KEYWORDS:
;       shift: if zero, behave normally by comparing the two specified
;       input spectra.  If non-zero, then ignore file2 and create a
;       fake second spectrum by simply shifting the wavelength scale
;       of the first spectrum.
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2011-Jan-06	GDW	Original version
;-

;; read in files...
data1 = mrdfits( file1, 'metadata', header)
spec1 = mrdfits( file1, 'spectrum', header)

data2 = mrdfits( file2, 'metadata', header)
spec2 = mrdfits( file2, 'spectrum', header)

;; determine whether to fake this out...
if ~ keyword_set(SHIFT) then shift=0.
if shift ne 0. then begin
    spec2.wavelength = spec1.wavelength + shift
    spec2.counts = spec1.counts
endif 

sep = ', '
print, data1.date, sep, $
       data1.grating, sep, $
       data1.gratepos, sep, $
       data1.central_wave, sep, $
       data1.blocking, sep, $
       data1.std_name
print, data2.date, sep, $
       data2.grating, sep, $
       data2.gratepos, sep, $
       data2.central_wave, sep, $
       data2.blocking, sep, $
       data2.std_name

;; overplot spectra...
offset = -1000.
plot, spec1.wavelength, spec1.counts, xtitle='Wavelength', ytitle='Counts', /nodata
oplot, spec1.wavelength, spec1.counts, color=gdwcolor('red')
oplot, spec2.wavelength, spec2.counts+offset, color=gdwcolor('blue')

;; determine offset...
delta = get_wavelength_shift( spec1.wavelength, spec1.counts, $
  spec2.wavelength, spec2.counts, status=status, chisq=chisq)
print, 'shift=', shift
print, 'delta=', delta

;; shift the spectrum...
spec2.wavelength += delta

;; overplot results...
offset = 2*offset
oplot, spec2.wavelength, spec2.counts+offset, color=gdwcolor('green')

end

;-----------------------------------------------------------------------
pro test_wavelength_shift, shift=shift
;-----------------------------------------------------------------------
; Front-end for testing get_wavelength_shift.
;
; EXAMPLE:
;       - log into HQ machine under dmoseng account
;       - launch IDL via ~dmoseng/bin/do_deimos_throughput
;       - in IDL, run this to compare two spectra
;               test_wavelength_shift
;       - in IDL, run this to compare spectrum against shifted version
;         of itself:
;               test_wavelength_shift, shift=500.
;-----------------------------------------------------------------------

thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq '' then $
  message, 'DEIMOS_THRU_DIR is not defined -- abort!'
extract_dir = thru_dir + '/extract/'

file1 = extract_dir + '1200G/2007nov11_d1111_0055.fits'
file2 = extract_dir + '1200G/2010jun14_d0614_0107.fits'

if ~ keyword_set(SHIFT) then shift=0.

doit, file1, file2, shift=shift
end

