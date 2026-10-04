# Generated PypeIt doc files that build-docs would change (MOIRCS)

build-docs was not run: it rewrites shared files in PypeIt/doc.  The
generators that depend on the spectrograph classes were run into a scratch
directory instead and compared with PypeIt/doc (branch
u/monodera/subaru_moircs_vbk at c08189b76):

```console
python regen_doc_tables.py <scratch dir>
```

1. doc/pypeit_par.rst (build_par_rst.py), MOIRCS section, 1 line
   (~l.11278):

   ```diff
   -          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits
   +          telgridfile = TellPCA_3000_26000_R10000.fits
   ```

   From the PCA telluric default (Implementation #7).  Nothing else
   changes: the VB_K values are in config_specific_par, which this
   listing does not show.

2. doc/include/inst_detector_table.rst (build_detector_table.py), the two
   subaru_moircs rows: the read noise is printed as 5.533985905294664
   instead of 5.534.  Same value (17.5/sqrt(10), from get_detector_par
   with hdu=None); only the formatting changes, because the value is now
   computed instead of hard-coded.  Possible cosmetic fix for prompt 10, in
   subaru_moircs.py: round the read noise (e.g. to 3 decimals) so the
   table keeps 5.534.

3. doc/include/spectrographs_table.rst (build_spectbl_rst.py): a
   subaru_moircs row would be added.  It is missing now (the HK500 work did
   not regenerate this table).  The file is also out of date on develop
   for other instruments: shane_hamspec is missing, the int_ids_eev10 row
   is joined to the end of the gtc_osiris_plus row, and the p200_tspec row
   differs.  So a rebuild changes ~175 lines, mostly column widths.

Unchanged: doc/include/data_dir.rst (build_cache_data_tbl.py lists
directories only; the new reid_arxiv file does not appear).

Other generators were not run (datamodels, bitmasks, script help,
dependencies, standards, telluric table): no MOIRCS change touches their
inputs.  No new script or datamodel.

Hand-written pages to replace with the drafts in this directory:
doc/spectrographs/subaru_moircs.rst and doc/tutorials/moircs_howto.rst.
Release notes: doc/releases/2.1.0dev.rst (see changelog.txt).
