;-----------------------------------------------------------------------
pro summarize_calib_files
;-----------------------------------------------------------------------
; Purpose:
;       Create file summary.txt which summarizes all of the files in
;       the current directory.
;-----------------------------------------------------------------------

;; get file list...
files = FILE_SEARCH('sens*.fits*')
out = 'summary.txt'
openw, ounit, out, /get_lun
;; format='(a34,x,a8,x,i5,x,a8,x,a10)'
;; format2='(a34,x,a8,x,a5,x,a8,x,a10)'
format  = '(%"%-34s %-8s %5d %-8s %-10s")'
format2 = '(%"%-34s %-8s %-5s %-8s %-10s")'
printf, ounit, '# filename', 'gratenam', 'wavelen', 'dwfilnam', 'targnam', format=format2

;; loop over files
n = n_elements(files)
for i=0,n-1 do begin

    ;; read header...
    f = files[i]
    h = mrdfits( f, 1, /silent)

    ;; print results...
    print, f, h.grating, nint(h.central_wave), h.blocking, h.std_name, format=format
    printf, ounit, f, h.grating, nint(h.central_wave), h.blocking, h.std_name, format=format
endfor 

free_lun, ounit
end 

