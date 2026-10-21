#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Append confidence assessment annotations (tier, tech/caller/tissue support, etc) as INFO fields for each variant
# ─────────────────────────────────────────────────────────────────────────────

. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod
module load bcftools/1.21 htslib/1.21

WORKDIR=$1
MERGE_MANIFEST="${WORKDIR}/merge_manifests/p25.merge.txt"
HEADER="${WORKDIR}/scripts/resources/header_for_collapse.txt"
COLUMNS="${WORKDIR}/scripts/resources/vcf_columns.txt"

awk -F',' 'NR > 1 {print $1,$2,$5}' ${MERGE_MANIFEST} | while read -r DONOR TISSUE_CODE INPUT_VCF
do 
    OUTPUT_VCF="${WORKDIR}/output/merged/annotated_vcfs/${DONOR}-${TISSUE_CODE}.vcf.gz"
    ANNOTATIONS="${WORKDIR}/output/merged/annotations/${DONOR}-${TISSUE_CODE}.tsv"

    bgzip -c ${ANNOTATIONS} > ${ANNOTATIONS}.gz
    tabix -f -s 1 -b 2 -e 2 ${ANNOTATIONS}.gz
    
    bcftools annotate --no-version ${INPUT_VCF} \
        -x INFO/AC \
        --annotations ${ANNOTATIONS}.gz \
        --header-lines ${HEADER} \
        --columns-file ${COLUMNS} \
        --threads 4 \
        --force \
    | bcftools view --no-version -Oz -o ${OUTPUT_VCF} 

    tabix -f -p vcf ${OUTPUT_VCF} 

done




