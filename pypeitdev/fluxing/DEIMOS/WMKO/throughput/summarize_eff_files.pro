;-----------------------------------------------------------------------
pro read_deimos_eff_sav_file, infile, extract, fextract, sensstd, efficiency
;-----------------------------------------------------------------------
extract = 0B
restore, infile
if size( extract, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named EXTRACT in file '+file
if size( fextract, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named FEXTRACT in file '+file
if size( sensstd, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named SENSSTD in file '+file
if size( efficiency, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named EFFICIENCY in file '+file

end

;-----------------------------------------------------------------------
pro summarize_eff_files
;-----------------------------------------------------------------------
; Purpose:
;       Create file summary.txt which summarizes all of the files in
;       the current directory.
;-----------------------------------------------------------------------

thru_dir = getenv('DEIMOS_THRU_DIR')
indir = thru_dir + '/eff/600ZD'
infiles = FILE_SEARCH(indir+'/*.sav')
ngood = intarr(4)
nbad = intarr(4)
passbands = intarr(4)

outfile = 'summarize_eff_files.dat'
openw, ounit, outfile, /get_lun

;; count files...
n = n_elements(infiles)
for i=0,n-1 do begin

    breakname, infiles[i], dirname, rootname, extn

    ;; read data....
    read_deimos_eff_sav_file, infiles[i], extract, fextract, $
      sensstd, efficiency

    date = extract.meta.date
    graname = extract.meta.grating
    wavelen = nint(extract.meta.central_wave) 
    targname = extract.meta.std_name

    ;; loop over available data...
    n_passbands = n_elements(efficiency)
    for j=0,n_passbands-1 do begin
        passband = nint(efficiency[j].lambda_eff)
        eff = efficiency[j].efficiency
        printf, ounit, rootname, graname, wavelen, targname, passband, eff, $
               format='(a,a6,i5,a12,i5,g12.5)'

        passbands[j] = passband
        if finite(eff) then ngood[j] += 1 else nbad[j] += 1
    endfor 

endfor 

for j=0,n_passbands-1 do $
  printf, ounit, passbands[j], ngood[j], nbad[j], $
         format='("Passband=",i, " Ngood=", i, " Nbad=", i)'

free_lun, ounit
end
