# User-defined fluxing parameters
[rdx]
   spectrograph = keck_nires
[fluxcalib]
   std_file = spec1d_s180604_0111-HIP107535_NIRES_2018Jun04T150456.336.fits
   sensfunc = HIP107535_keck_nires.fits
   # The next 2 are make-believe
   star_type = A0
   star_mag = 10.


flux read
  spec1d_s180604_0089-J1724+1901_NIRES_2018Jun04T125612.386.fits J1217+3905_1.fits
flux end