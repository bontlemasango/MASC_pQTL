// Necessary folders and commands
conda='eval "\$(/home/b/bl212/y/micromamba shell hook --shell bash)"'
activate="${conda}; micromamba activate regenie"
data_dir="../../data"
raw="${data_dir}/raw"
intermediate="${data_dir}/intermediate"
processed="${data_dir}/processed"

process get_proteins{
    executor="local"
    publishDir "$intermediate", mode:"rellink"
    output:
        path "proteins.tsv", emit: "pheno"
        path "proteomics_map.tsv", emit: "map"
    script:
    """
    ${activate}
    get_proteins.R \
        --output proteins.tsv \
        --mapping-file proteomics_map.tsv
    """
}

process get_covars{
    executor="local"
    publishDir "$intermediate", mode:"rellink"
    input:
        path proteins
        path map
    output:
        path "proteins.covars.tsv", emit: "covars"
    script:
    """
    ${activate}
    get_covars.R ${proteins} ${map} \
        --output proteins.covars.tsv
    """
}

process get_GWAS_covars{
    executor="local"
    publishDir "$intermediate", mode:"rellink"
    output:
        path "exceed.eur.pcs.tsv"
    script:
    """
    ${activate}
    get_GWAS_covar.R \
        --output exceed.eur.pcs.tsv
    """
}

process residualise{
    executor="local"
    publishDir "$intermediate", mode:"rellink"
    input:
        path proteins
        path covars
    output:
        path "proteins.resid.tsv", emit: "pheno"
    script:
    """
    ${activate}
    residualise.R ${proteins} ${covars} \
        --output proteins.resid.tsv
    """
}

process inverse_normal_transform{
    executor="local"
    publishDir "$intermediate", mode:"rellink"
    input:
        path proteins
    output:
        path "proteins.resid.INT.tsv", emit: "pheno"
    script:
    """
    ${activate}
    rank_inverse_normal.R ${proteins} \
        --output proteins.resid.INT.tsv
    """
}


process create_regenie_script{
    executor="local"
    publishDir "$processed", mode:"copy"
    input:
        path template_file1
        path template_file2
        val proteins
        val covars
    output:
        path "STEP_1-RegenieStep1.sbatch"
        path "STEP_2-RegenieStep2.sbatch"
    script:
    """
    ${activate}
    inject.py ${template_file1} \
        --sub JOB_NAME=pQTL_regenie_step1 \
              PHENO=${proteins} \
              COVAR=${covars} \
              OUTPUT=regenie_step1 \
        --output STEP_1-RegenieStep1.sbatch
    
    inject.py ${template_file2} \
        --sub JOB_NAME=pQTL_regenie_step2 \
              PHENO=${proteins} \
              COVAR=${covars} \
              STEP1_OUT=regenie_step1 \
        --output STEP_2-RegenieStep2.sbatch
    """
}


workflow {
    proteins = get_proteins()
    covars = get_covars(proteins.pheno, proteins.map)
    geno_covars = get_GWAS_covars()
    
    proteins = residualise(proteins.pheno, covars.covars) \
                | inverse_normal_transform
    
    create_regenie_script(file("${raw}/regenie_step1.slurm.inject"),
                          file("${raw}/regenie_step2.slurm.inject"),
                          proteins.pheno, geno_covars)
}
