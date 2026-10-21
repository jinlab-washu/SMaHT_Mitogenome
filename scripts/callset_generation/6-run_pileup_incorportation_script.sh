#!/bin/bash
# ────────────────────────────────────────────────────────────────────────────
# Generate updated donor-tissue VCFs with pileup info incorporated. 
# Pileup VCFs are generated for each donor-tissue individually for each sequencing
# technology using minipileup. Individual pileup VCFs are merged across samples  
# for each technology and the three merged VCFs (illumina, ont, pacbio) are provided
# as input to the merge_caller_and_pileup.py script alongside each original 
# donor-tissue VCF.
# ─────────────────────────────────────────────────────────────────────────────

module load modules modules-init modules-gs
module load bcftools/1.21

INDIR=$1
PILEUP_PACBIO=$2
PILEUP_ONT=$3
PILEUP_ILLUMINA=$4

VCFS=$(ls ${INDIR}/*.vcf.gz)

mkdir -p ${INDIR}/pileup/

for VCF in ${VCFS[@]};
do
    sample=$(basename ${VCF} .vcf.gz)
    echo $sample
    
    merge_caller_and_pileup.py \
        -s ${sample} \
        -c ${INDIR}/${sample}.vcf.gz \
        -o ${INDIR}/pileup/${sample}.withpileup.vcf.gz \
        -p ${PILEUP_PACBIO} \
        -n ${PILEUP_ONT} \
        -i ${PILEUP_ILLUMINA}

    bcftools annotate -x INFO/TechSupport,INFO/ToolSupport,INFO/Illumina_AF,INFO/PacBio_AF,INFO/ONT_AF \
        ${INDIR}/pileup/${sample}.withpileup.vcf.gz -Oz -o ${INDIR}/pileup/${sample}.withpileup.final.vcf.gz
done

