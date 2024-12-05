# pQTL mapping using REGENIE

## Programs used
This script uses three pieces of software: `nextflow`, `R` and `python`.  
Nextflow can be installed using the following command:  
`curl -s https://get.nextflow.io | bash`  
Further information on Nextflow and it's installation can be found [here](https://www.nextflow.io/docs/latest/install.html).

To install the R and python packages used, please use the `environment.yml` file in `code/`, eg:
`conda env create -f environment.yml`

## The workflow
Due to some limitations of the servers used, background tasks like nextflow cannot be run for longer than 2 hours. Because of this, the nextflow script runs the phenotype preparation and creates regenie scripts with the appropriate inputs to be used by the `SLURM Workload Manager` and run the proteomic GWAS.

If the phenotype files are already prepared, it may be that only the final step in the nextflow workflow - `create_regenie_script` is needed. Use the files in `./data/raw` with extension `.slurm.inject` as templates for your regenie scripts in this case.

If using the nextflow pipeline, most of the scripts used are in `./code/proteomic_GWAS/bin/`. All of these files should be modified using `chmod +x FILENAME_HERE` to be executable by nextflow.

The final command at the bottom of the nextflow script, `workflow{...}`, is the main script that is run. If steps should be skipped, please alter this script by removing the appropriate steps you would not like to be run.

Changes from the original:
  - Files for the ALICE filesystem must be copied from /rfs to /scratch before running the regenie script, and subsequently deleted
  - PCA creation and analysis needs to be added (use plink/genotyped)
  - Some files will need to be mapped to their appropriate barcode values \(\_b0 or \_b123\)