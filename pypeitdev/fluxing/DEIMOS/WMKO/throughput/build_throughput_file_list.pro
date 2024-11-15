;-----------------------------------------------------------------------
pro build_throughput_file_list
;-----------------------------------------------------------------------
;+
; NAME:
;	BUILD_THROUGHPUT_FILE_LIST
;
; PURPOSE:
;	This routine will extract critical header parameters and put
;	them into a catalog file with info on 
;       filename, grating, cenlam, filter, and star
;
; CATEGORY:
;	IPM
;
; CALLING SEQUENCE:
;	build_throughput_file_list
;
; INPUTS:
;       None
;
; OUTPUTS:
;	Creates output file 'input.lst' which is used as input for
;	run_review_deimos_throughput
;
; RESTRICTIONS:
;
; EXAMPLE:
;       1) To generate list, first launch IDL with the correct
;       setup by executing this in the shell:
;               do_deimos_throughput
;       then at the IDL prompt type
;               build_throughput_file_list
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2011-Oct-31	GDW	Original version
;       2017-Dec-15     CA      Create a copy of the file input.lst in
;                               $DEIMOS_THRU_DIR/dat/fiduclias.csv
;                               Fixed bug with informational messages.
;       2017-Dec-28     CA      Define two different structures, one
;                               for $DEIMOS_THRU_DIR/dat/input.lst and
;                               the other one for the file
;                               $DEIMOS_THRU_DIR/dat/fiducials.csv
;                               The structure for 
;                               fiducial_file.pro
;-
;-----------------------------------------------------------------------

;; define files...
thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq '' then $
  message, 'DEIMOS_THRU_DIR is not defined -- abort!'


;; list of approved standard star names in
;; $XIDL_DIR/Spec/Flux/x_stdcalibfil.pro that
;; will be used in the script. No matter what
;; target name format is used in the 
;; observe_flux_standard script, the target name
;; that appears in the file
;; $DEIMOS_THRU_DIR/dat/fiducials.csv should match
;; one of the names in the following list.

standards = [ 'Feige34', 'Feige110', $
              'BD+28deg4211', 'BD+33deg2642', $
              'HZ_44', 'G191B2B']

data_dir = thru_dir + '/dat/'
raw_dir = thru_dir + '/raw/'
all_list = data_dir + 'all.lst'
input_list = data_dir + 'input.lst'
fiducials_csv = data_dir + 'fiducials.csv'

;; get list of files...
;; readcol, all_list, format='(a)', files
q = '''
;; '

command = 'find ' + thru_dir + '/raw -name ' + $
          q + '*.fits*' + q + ' -print | fgrep -v BAD'
print, command
spawn, command, files
n_files = n_elements(files)

;; allocate output struct for input.lst
s = { filename:'', targname:'', gratenam:'', wavelen:0, dwfilnam:'', $
     dateobs:'', status:1B }
info = replicate(s, n_files)

;; allocate output struct for fiducails.csv
s_fid = { gratenam:'', dwfilnam:'', wavelen:0, targname:'', dir:'', $
     filenam:'', dateobs:'',  status:1B }
info_fid = replicate(s_fid, n_files)

;; scan through files...
for i=0,n_files-1 do begin

   ;; default status is "good"...
   t = s
   t_fid = s_fid
   
   ;; read file header...
   file = files[i]
   print, '[', + stringify(i) + '/' + stringify(n_files) + '] file=', file
   hdr = headfits( file, exten=0)
   
   ;; check keywords...
   t.filename = file
   t.dateobs  = strtrim(sxpar( hdr, 'DATE-OBS'), 2)
   t.targname = strtrim(sxpar( hdr, 'TARGNAME'), 2)
   t.dwfilnam = strtrim(sxpar( hdr, 'DWFILNAM'), 2)
   t.gratenam = strtrim(sxpar( hdr, 'GRATENAM'), 2)
   gratepos = sxpar( hdr, 'GRATEPOS')
   if gratepos eq 3 then begin
      keyword = 'G3TLTWAV'
      t.wavelen = nint(sxpar( hdr, keyword))
   endif else if gratepos eq 4 then begin
      keyword = 'G4TLTWAV'
      t.wavelen = nint(sxpar( hdr, keyword))
   endif else begin
      message, 'WARNING: invalid GRATEPOS (' + $
               strtrim(string(gratepos), 2) + $
               ') in file '+file, /info
      t.status = 0B
   endelse 

   ;; need to split the filename into file and path 
   ;; for the fiducials_csv struct.
   
   filename_split = strsplit(file, '/', /EXTRACT)
   n_split = n_elements(filename_split)
   filename_only = filename_split[n_split-1]
   filename_only_split = strsplit(filename_only, '.', /EXTRACT)
   filename_root = filename_only_split[0]
   file_path = filename_split[n_split-2]

   ;; find the XIDL standard star name that matches the
   ;; image header target name
   
   if stregex(t.targname, '34', /BOOLEAN) then xstandard = standards[0]
   if stregex(t.targname, '110', /BOOLEAN) then xstandard = standards[1]
   if stregex(t.targname, '4211', /BOOLEAN) then xstandard = standards[2]
   if stregex(t.targname, '2642', /BOOLEAN) then xstandard = standards[3]
   if stregex(t.targname, '44', /BOOLEAN) then xstandard = standards[4]
   if stregex(t.targname, '191', /BOOLEAN) then xstandard = standards[5]

   t_fid.gratenam = t.gratenam
   t_fid.dwfilnam = t.dwfilnam
   t_fid.wavelen = t.wavelen
   t_fid.targname = xstandard
   t_fid.dir = file_path
   t_fid.filenam = filename_root
   t_fid.dateobs = t.dateobs

   ;; flag invalid dates...
   if stregex( t.dateobs, $
               '^[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]$', $
               /boolean) ne 1 then begin
      message, 'WARNING: invalid date in file '+files[i], /info
      t.status = 0B
   endif

   t_fid.status = t.status
   
   ;; insert into struct for input.lst
   info[i] = t

   ;; insert into struct for fiducials.csv
   info_fid[i] = t_fid
   
endfor 

;; remove invalid dates from the input.lst struct
good = where( info.status eq 1B )
info = info[good]

;; sort by date the input.lst struct
order = sort(info.dateobs)
info = info[order]

;; remove invalid dates from the input.lst struct
good_fid = where( info_fid.status eq 1B )
info_fid = info_fid[good_fid]

;; sort by date the input.lst struct
order_fid = sort(info_fid.dateobs)
info_fid = info_fid[order_fid]

;; write to file input.lst...

write_csv, input_list, info, $
           header=['filename', 'targname', 'gratenam', 'wavelen', $
                   'dwfilnam', 'dateobs', 'status' ]

message, 'Wrote file '+input_list, /info

;; write to file fiducials.csv...

write_csv, fiducials_csv, info_fid, $
           header=['gratenam', 'dwfilnam', 'wavelen', 'targets', $
                   'dir', 'filename', 'dateobs', 'status']

message, 'Wrote file '+fiducials_csv, /info

end

