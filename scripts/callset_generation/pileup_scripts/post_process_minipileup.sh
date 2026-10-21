#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Pileup VCFs are processed using bcftools to normalize multiallelics, add 
# additional format fields (DPF,DPR,AF), and separate snvs and indels. Each
# individual pileup VCF is then processed using annotate_pileup.py for 
# strand bias and coverage-based alt read support checks for each variant call.
# Final pileup VCFs are merged across samples and snvs + indels VCFs are
# concatenated.
# ─────────────────────────────────────────────────────────────────────────────

module load modules modules-init modules-gs
module load bcftools/1.21 htslib/1.21

### set these variables
TECH=$1
OUTDIR=$2
###
SCRIPTDIR="/net/nwgc/vol1/home/czaka/analysis/minipileup/scripts/"
chrom_map="${SCRIPTDIR}/chrom_map.txt"
REF="/net/nwgc/vol1/home/czaka/ref/MT.fasta"
MAX_CONCURRENT_JOBS=24

mkdir -p ${OUTDIR}/${TECH}

for VCF in $(ls ${OUTDIR}/${TECH}/*.minipileup.vcf);
do  
    SAMPLE_NAME=$(basename ${VCF} .minipileup.vcf)
    echo $SAMPLE_NAME
    
    OLD_NAME=$(bcftools query -l "$VCF")

    bcftools +fill-tags ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.vcf -- -t FMT/DPF:1='smpl_sum(FMT/ADF)',FMT/DPR:1='smpl_sum(FMT/ADR)' \
    | bcftools annotate -x QUAL --rename-chrs ${chrom_map} \
    | bcftools norm -f ${REF} --multiallelics -any \
    | bcftools norm -f ${REF} --atomize --atom-overlaps . \
    | bcftools view -i 'FMT/ADF[0:1]>=2 & FMT/ADR[0:1]>=2' \
    | bcftools +fill-tags  -- -t FMT/AF:1='(FMT/ADF[:1]+FMT/ADR[:1])/(FMT/DPF+FMT/DPR)' \
    | bcftools reheader -s <(echo "${OLD_NAME} ${SAMPLE_NAME}") \
    | bcftools view -o ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.vcf

    bcftools view -v snps ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.vcf -Oz -o ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.snps.vcf.gz
    bcftools view -v indels ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.vcf -Oz -o ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.indels.vcf.gz

    tabix -p vcf ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.snps.vcf.gz
    tabix -p vcf ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.indels.vcf.gz

    ${SCRIPTDIR}/annotate_pileup.py \
        -i ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.snps.vcf.gz \
        -o ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.snps.anno.vcf.gz > ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.snps.anno.log 2>&1 &

    ${SCRIPTDIR}/annotate_pileup.py \
        -i ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.indels.vcf.gz \
        -o ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.indels.anno.e01.vcf.gz -e 0.01 > ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.indels.anno.e01.log 2>&1 &

    # allow to execute up to max n jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge ${MAX_CONCURRENT_JOBS} ]]; then
        wait -n
    fi

done

wait

for VCF in $(ls ${OUTDIR}/${TECH}/*.minipileup.vcf);
do  
    SAMPLE_NAME=$(basename ${VCF} .minipileup.vcf)
    echo $SAMPLE_NAME

    tabix -f -p vcf ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.snps.anno.vcf.gz
    tabix -f -p vcf ${OUTDIR}/${TECH}/${SAMPLE_NAME}.minipileup.norm.indels.anno.e01.vcf.gz

done

LIST_OF_VCFS=$(ls ${OUTDIR}/${TECH}/*.minipileup.norm.snps.anno.vcf.gz)
bcftools merge -m none -M ADF:0,ADR:0 --force-single ${LIST_OF_VCFS} -o ${OUTDIR}/minipileup.merged.${TECH}.snps.vcf

LIST_OF_VCFS=$(ls ${OUTDIR}/${TECH}/*.minipileup.norm.indels.anno.e01.vcf.gz)
bcftools merge -m none -M ADF:0,ADR:0 --force-single ${LIST_OF_VCFS} -o ${OUTDIR}/minipileup.merged.${TECH}.indels.e01.vcf

bgzip -c ${OUTDIR}/minipileup.merged.${TECH}.snps.vcf > ${OUTDIR}/minipileup.merged.${TECH}.snps.vcf.gz
bgzip -c ${OUTDIR}/minipileup.merged.${TECH}.indels.e01.vcf > ${OUTDIR}/minipileup.merged.${TECH}.indels.e01.vcf.gz

tabix -f -p vcf ${OUTDIR}/minipileup.merged.${TECH}.snps.vcf.gz
tabix -f -p vcf ${OUTDIR}/minipileup.merged.${TECH}.indels.e01.vcf.gz

bcftools concat --allow-overlaps ${OUTDIR}/minipileup.merged.${TECH}.snps.vcf.gz ${OUTDIR}/minipileup.merged.${TECH}.indels.e01.vcf.gz \
    -Oz -o ${OUTDIR}/minipileup.merged.${TECH}.vcf.gz

tabix -f -p vcf ${OUTDIR}/minipileup.merged.${TECH}.vcf.gz




