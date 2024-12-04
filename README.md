# pQTL mapping using REGENIE

Changes from the original:
  - Files for the ALICE filesystem must be copied from /rfs to /scratch before running the regenie script, and subsequently deleted
  - PCA creation and analysis needs to be added (use plink/genotyped)
  - Some files will need to be mapped to their appropriate barcode values \(\_b0 or \_b123\)