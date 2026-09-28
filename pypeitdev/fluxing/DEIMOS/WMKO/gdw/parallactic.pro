;-----------------------------------------------------------------------
function parallactic, ha, dec, lat
;-----------------------------------------------------------------------
; Compute the parallactic angle as a function of hour angle (HA) and
; declination (Dec) for a given observer latitude (lat).  Note
; that HA is expected in hours ( <0 for East and >0 for West) and Dec
; and Lat are expected in degrees.  Parallactic is returned in deg.
;-----------------------------------------------------------------------

;; define vector from observer to target...
ha_deg = -ha*15.d0
T = RD2XYZ( ha_deg, dec)

;; define vector from observer to zenith...
ha_deg = 0.d0
Z = RD2XYZ( ha_deg, lat)

;; define vector (B) from target to zenith...
A = cross(Z,T)
B = cross(T,A)

;; define vector (N) to north...
ha_north = 0.d0
dec_north = 90.d0
N = RD2XYZ( ha_north, dec_north)

;; define vector (D) from target to north...
C = cross(N,T)
D = cross(T,C)

;; parallactic angle is defined as the position angle of B w.r.t. D at
;; T...
parang = PA( B, D, T)

;; range checking...
if dec gt lat then begin
    parang_max = 360.
endif else begin
    parang_max = 180.
endelse

parang = parang + 720.
while parang gt parang_max do begin
    parang = parang - 360.
endwhile

return, parang

end
