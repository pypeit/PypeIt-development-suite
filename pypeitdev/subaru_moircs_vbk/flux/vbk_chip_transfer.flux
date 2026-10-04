# Option B: one spec1d per detector, each with its own sensitivity function.
# DET01: the chip-1 sensfunc (HIP59174 std1 + std2, A0 Vega, PCA tellurics).
# DET02: two copies of the std1 A DET02 objects (the science pair has no
# DET02 object), one per transfer mode.
flux read
  path /mnt/s1data01/work/research/pypeit-development/PypeIt-development-suite/pypeitdev/subaru_moircs_vbk/flux/fluxed
  filename | sensfile
  spec1d_MCSA00352031-CC0958_PA200_1_MOIRCS_20260529T060511.128_DET01.fits | /mnt/s1data01/work/research/pypeit-development/PypeIt-development-suite/pypeitdev/subaru_moircs_vbk/flux/sens_A0_pca/sens_HIP59174_A0_pca.fits
  spec1d_MCSA00352033-CC0958_PA200_1_MOIRCS_20260529T060843.130_DET01.fits | /mnt/s1data01/work/research/pypeit-development/PypeIt-development-suite/pypeitdev/subaru_moircs_vbk/flux/sens_A0_pca/sens_HIP59174_A0_pca.fits
  spec1d_MCSA00352125-HIP59174_MOIRCS_20260529T082520.088_DET01.fits | /mnt/s1data01/work/research/pypeit-development/PypeIt-development-suite/pypeitdev/subaru_moircs_vbk/flux/sens_A0_pca/sens_HIP59174_A0_pca.fits
  spec1d_MCSA00352125-HIP59174_MOIRCS_20260529T082520.088_DET02_direct.fits | /mnt/s1data01/work/research/pypeit-development/PypeIt-development-suite/pypeitdev/subaru_moircs_vbk/flux/sens_HIP59174_A0_pca_det2_direct.fits
  spec1d_MCSA00352125-HIP59174_MOIRCS_20260529T082520.088_DET02_flatratio.fits | /mnt/s1data01/work/research/pypeit-development/PypeIt-development-suite/pypeitdev/subaru_moircs_vbk/flux/sens_HIP59174_A0_pca_det2_flatratio.fits
flux end
