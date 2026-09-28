;-----------------------------------------------------------------------
function parse_review_database, db, infile
;-----------------------------------------------------------------------
; purpose: return review database record matching specified key field 
; (infile)
;-----------------------------------------------------------------------

;; clean the filename...
in = trim_gz(strtrim(infile,2))

;; locate appropriate record...
n = n_elements(db)
for i=0,n-1 do begin
    if strtrim(db[i].filename,2) eq in then begin

        ;; trim leading and trailing spaces...
        record = db[i]
        record.filename = strtrim(record.filename,2)
        record.status   = strtrim(record.status,2)
        record.comment  = strtrim(record.comment,2)
        return, record
    endif 
endfor 

;; check for no match...
message, 'WARNING: no match in review database for filename '+in, /info
return, -1

end

