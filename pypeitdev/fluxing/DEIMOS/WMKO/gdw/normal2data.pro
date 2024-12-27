;-----------------------------------------------------------------------
pro normal2data, x_norm, y_norm, x_data, y_data
;-----------------------------------------------------------------------
; Purpose:
;	Returns the data coordinates corresponding the normalized
;	coordinates for the current plot.
;-----------------------------------------------------------------------

x_min = !x.crange[0]
x_max = !x.crange[1]
y_min = !y.crange[0]
y_max = !y.crange[1]

x_data = x_min + x_norm*(x_max-x_min)
y_data = y_min + y_norm*(y_max-y_min)

end


