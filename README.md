# pQTL mapping using REGENIE

## Programs used
This script uses three pieces of software: `nextflow`, `R` and `python`.  
Nextflow can be installed using the following command:  
`curl -s https://get.nextflow.io | bash`  
Further information on Nextflow and it's installation can be found [here](https://www.nextflow.io/docs/latest/install.html).

## Installation of conda

The script uses micromamba, a faster version of conda. This can be isntalled in a Linux environment using the following command: `"${SHELL}" <(curl -L micro.mamba.pm/install.sh)`  
Detailed installation instructions here: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

For the more standard approach install miniconda instead. If you use this one, replace any code which says 'micromamba' with 'conda'. Instructions to install here: https://docs.anaconda.com/miniconda/install/#quick-command-line-install  
Don't install both, you only need one!

## Regenie/R/Python setup
To install all the necessary components of the pipeline (including REGENIE) use the following:
1. cd into the GitHub folder proteomic_GWAS:  
1. `cd code`
1. `micromamba env create -f environment.yml`

Everytime you want to use the workflow, use `micromamba activate regenie`. You can check this worked by seeing if 'regenie' autocompletes in your command line.

## The workflow
Due to some limitations of the servers used, background tasks like nextflow cannot be run for longer than 2 hours. Because of this, the nextflow script runs the phenotype preparation and creates regenie scripts with the appropriate inputs to be used by the `SLURM Workload Manager` and run the proteomic GWAS.

If the phenotype files are already prepared, it may be that only the final step in the nextflow workflow - `create_regenie_script` is needed. Use the files in `./data/raw` with extension `.slurm.inject` as templates for your regenie scripts in this case.

If using the nextflow pipeline, most of the scripts used are in `./code/proteomic_GWAS/bin/`. All of these files should be modified using `chmod +x FILENAME_HERE` to be executable by nextflow.

The final command at the bottom of the nextflow script, `workflow{...}`, is the main script that is run. If steps should be skipped, please alter this script by removing the appropriate steps you would not like to be run.

Changes from the original:
  - Files for the ALICE filesystem must be copied from /rfs to /scratch before running the regenie script, and subsequently deleted
  - PCA creation and analysis needs to be added (use plink/genotyped)
  - Some files will need to be mapped to their appropriate barcode values \(\_b0 or \_b123\)