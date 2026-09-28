pro make_review_database
;+
; NAME:
;	MAKE_REVIEW_DATABASE
;
; PURPOSE:
;	This program will merge the indepdent files good.lst, bad.lst,
;	and review.lst into a single file review.dat
;
; CATEGORY:
;	IPM
;
; CALLING SEQUENCE:
;	make_review_database
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2011-Nov-15	GDW	Original version
;-

verbose = 0B

;; define files...
thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq '' then $
  message, 'DEIMOS_THRU_DIR is not defined -- abort!'
data_dir = thru_dir + '/dat/'
input_list   = data_dir + 'input.lst'
reject_list = data_dir + 'bad.lst'
approve_list = data_dir + 'good.lst'
review_list = data_dir + 'review.lst'
output_list = data_dir + 'review.fits'

;; get input list...
result = read_csv( input_list)
all_files = result.field1
n_files = n_elements(all_files)

;; get approved list...
readcol, approve_list, good_files, format='(a)'
n_good = n_elements(good_files)

;; get bad list...
result = read_csv(reject_list)
bad_files = result.field1
bad_comment = result.field2
n_bad = n_elements(bad_files)

;; get review list...
result = read_csv(review_list)
review_files = result.field1
review_comment = result.field2
n_review = n_elements(review_files)

;; fix filenames..
all_files = trim_gz(all_files)
bad_files = trim_gz(bad_files)
good_files = trim_gz(good_files)
review_files = trim_gz(review_files)

;; generate a struct...
a = {filename:'', status:'', timestamp:'', comment:''}
data = replicate( a, n_files)

;; insert data into file...
for i=0,n_files-1 do begin

    filename = all_files[i]

    ;; insert filename...
    data[i].filename = filename

endfor 

;; check for good...
for j=0,n_good-1 do begin
    found = 0B

    for i=0,n_files-1 do begin
        if good_files[j] eq all_files[i] then begin
            found = 1B

            if data[i].status ne '' then begin
                print, 'WARNING: good file ', all_files[i], $
                       ' already marked as ', data[i].status
            endif else begin
                data[i].status = 'good'
                if verbose then print, 'marked good on ', all_files[i]
            endelse 
        endif
    endfor 

    if ~ found then message, 'ERROR: no match for good file '+good_files[j]

endfor 

;; check for bad...
for j=0,n_bad-1 do begin
    found = 0B
    
    for i=0,n_files-1 do begin
        if bad_files[j] eq all_files[i] then begin
            found = 1B
            if data[i].status ne '' then begin
                print, 'WARNING: bad file ', all_files[i], $
                       ' already marked as ', data[i].status
            endif else begin
                data[i].status = 'bad'
                data[i].comment = bad_comment[j]
                if verbose then print, 'marked bad on ', all_files[i]
            endelse 
        endif
    endfor 

    if ~ found then message, 'ERROR: no match for bad file '+bad_files[j]

endfor 

;; check for review...
for j=0,n_review-1 do begin
    found = 0B
    
    for i=0,n_files-1 do begin
        if review_files[j] eq all_files[i] then begin
            found = 1B
            if data[i].status ne '' then begin
                print, 'WARNING: review file ', all_files[i], $
                       ' already marked as ', data[i].status
            endif else begin
                data[i].status = 'review'
                data[i].comment = review_comment[j]
                if verbose then print, 'marked review on ', all_files[i]
            endelse 
        endif
    endfor 

    if ~ found then message, 'ERROR: no match for review file '+review_files[j]

endfor 

;; save file...
mwrfits, data, output_list, /create

end
