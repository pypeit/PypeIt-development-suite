;-----------------------------------------------------------------------
pro throughput_update_list, filename, v1, v2
;-----------------------------------------------------------------------

;; read data...
if n_params() eq 2 then begin

    ;; if file doesn't exist, then just create it...
    if ~ file_test(filename) then begin
        openw, ounit, filename, /get_lun
        printf, ounit, v1
        free_lun, ounit
    endif else begin

        ;; get current contents...
        readcol, filename, format='(a)', c1, /silent
        match = where( c1 eq v1, count)
        
        if count eq 0 then begin
            
            ;; add new value to file...
            openw, ounit, filename, /get_lun, /append
            printf, ounit, v1
            free_lun, ounit
            
        endif else begin
            message, 'WARNING: '+v1+' already appears in file '+filename, /info
        endelse 
    endelse 

endif else if n_params() eq 3 then begin

    ;; if file doesn't exist, then just create it...
    if ~ file_test(filename) then begin
        openw, ounit, filename, /get_lun
        printf, ounit, v1, v2
        free_lun, ounit
    endif else begin

        ;; get current contents...
        readcol, filename, format='(a,f)', c1, c2, /silent
        
        match = where( c1 eq v1, count)
        
        if count eq 0 then begin
            
            ;; add new values to arrays...
            c1 = [c1,v1]
            c2 = [c2,v2]
            
        endif else if count eq 1 then begin
            
            ;; update existing array value...
            c2[match] = v2
            
        endif else begin
            message, 'WARNING: '+v1+' appears multiple times in file '+filename, $
                     /info
            return
        endelse 
        
        openw, ounit, filename, /get_lun, width=132
        for i=0,n_elements(c1)-1 do printf, ounit, c1[i], c2[i]
        free_lun, ounit
    endelse 
        
endif else begin
    message, 'wrong number of parameters'
endelse 

end
