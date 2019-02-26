# User-defined fluxing parameters
[rdx]
   spectrograph = shane_kast_blue
[fluxcalib]
   std_file = spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits
   sensfunc = feige66_shane_kast_blue.fits

flux read
  spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits J1217+3905_1.fits
  spec1d_b28-J1217p3905_KASTb_2015May20T051801.470.fits J1217+3905_2.fits
flux end