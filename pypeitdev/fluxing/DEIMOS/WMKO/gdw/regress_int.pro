;-----------------------------------------------------------------------
pro regress_int, x, y, slope, yfit, intercept
;-----------------------------------------------------------------------
; regress_int.awk / GDW / 7 Dec 1991
;
; This routine computes the x-axis intercept of the line with the specified
; slope (slope) which minimizes the vertical deviation.  
;
;-----------------------------------------------------------------------

nx = n_elements(x)
ny = n_elements(y)
if nx ne ny then begin
    message, 'nx and ny differ'
endif

intercept = (total(y) - slope*total(x))/nx
yfit = slope*x + intercept

return
end
