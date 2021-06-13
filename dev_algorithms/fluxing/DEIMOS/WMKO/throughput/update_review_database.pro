;-----------------------------------------------------------------------
pro update_review_database, filename, db, infile, $
  STATUS=status, COMMENT=comment
;-----------------------------------------------------------------------
; Purpose: update a record in the review database
;
; Args:
;       filename: name of database file on disk [I]
;       db: struct with data [I]
;       infile: name of record to update (key) [I]
;-----------------------------------------------------------------------

;; clean the filename...
in = trim_gz(strtrim(infile,2))

;; locate appropriate record...
found = 0B
n = n_elements(db)
for i=0,n-1 do begin
    if strtrim(db[i].filename,2) eq in then begin
        found = 1B

        ;; update specified fields...
        if keyword_set(STATUS) then db[i].status = status
        if keyword_set(COMMENT) then db[i].comment = comment

        ;; update timestamp...
        db[i].timestamp = systime(0)
        break
    endif 
endfor 

;; check for no match...
if ~ found then begin
    message, 'WARNING: no match in review database for filename '+in, /info
    return
endif

;; update on disk...
mwrfits, db, filename, /create

end

