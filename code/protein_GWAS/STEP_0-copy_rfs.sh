# Script to copy the rfs files where necessary
TEMP_FOLDER=/scratch/gen1/bl212/exceed/proteomic_GWAS_freeze2/data/raw/tmp

PLINK_FILE=/rfs/EXCEED/Genotyping/freeze2/plink/genotyped/EXCEED_freeze2_genotyped_TOPMed_b38
BGEN_FILE=/rfs/EXCEED/Genotyping/freeze2/bgen/chr22.exceed.freeze2.topmed.gp.8bit.varids

# Copy the plink files to the temporary folder
for extension in .bed .bim .fam; do
    TEMP_FILE=${TEMP_FOLDER}/$(basename $PLINK_FILE)${extension}
    if [ ! -f $TEMP_FILE ]; then
        cp ${PLINK_FILE}${extension} $TEMP_FILE
        echo copied to $TEMP_FILE
    fi
done

# Repeat for bgen files
for i in {1..23}; do
    for extension in .bgen .bgen.bgi .bgen.sample; do
        FILE_NAME=$(basename $BGEN_FILE)
        FILE_NAME=${FILE_NAME//chr22/chr${i}}${extension}
        TEMP_FILE=${TEMP_FOLDER}/$FILE_NAME
        if [ ! -f ${TEMP_FILE} ]; then
            cp ${BGEN_FILE//chr22/chr${i}}${extension} ${TEMP_FILE}
            echo copied to ${TEMP_FILE}
        fi
    done
done