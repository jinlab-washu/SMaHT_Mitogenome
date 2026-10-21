#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Format input VCFs across callers (mutect2, mutserve, mitoscope, himito) to standardize prior to merging
# Includes standardizing sample names, chromosome names, genotypes, info/format fields
# Remove blacklisted sites (m.3107), normalize multiallelics, filter for PASS variants
# ─────────────────────────────────────────────────────────────────────────────

module load modules modules-init modules-gs
module load bcftools/1.21 htslib/1.21

MUTECT2_MANIFEST=$1
MUTSERVE_MANIFEST=$2
MITOSCOPE_PB_MANIFEST=$3
MITOSCOPE_ONT_MANIFEST=$4
HIMITO_PB_MANIFEST=$5
HIMITO_ONT_MANIFEST=$6

WORKDIR=$7
OUTDIR="${WORKDIR}/output"
MERGE_MANIFEST="${WORKDIR}/merge_manifests/p25.txt"
chrom_map="${WORKDIR}/scripts/resources/chrom_map.txt"
REF="${WORKDIR}/scripts/resources/MT.fasta"

echo "sample_id,file" > ${MERGE_MANIFEST}

## mutect2 illumina
awk -F',' '{print $1,$2}' ${MUTECT2_MANIFEST} | while read -r SAMPLE VCF_FILE
do
    tech='illumina'
    echo "mutect2"
    echo $SAMPLE
    mkdir -p ${OUTDIR}/${tech}

    OLD_NAME=$(bcftools query -l "$VCF_FILE")
    NEW_NAME="${SAMPLE}-mutect2"

    bcftools reheader -s <(echo "${OLD_NAME} ${NEW_NAME}") ${VCF_FILE} \
    | bcftools +setGT -- -t a -n c:'1' \
    | bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/AF \
    | bcftools filter -e 'FILTER != "PASS" && FMT/AF < 0.95' --mode x -Oz -o ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz

    bcftools index --tbi ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz
    echo "${NEW_NAME},${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz" >> ${MERGE_MANIFEST}
done

## mutserve illumina
awk -F',' '{print $1,$2}' ${MUTSERVE_MANIFEST} | while read -r SAMPLE VCF_FILE
do
    echo "mutserve"
    tech='illumina'
    echo $SAMPLE
    mkdir -p ${OUTDIR}/${tech}

    OLD_NAME=$(bcftools query -l "$VCF_FILE")
    NEW_NAME="${SAMPLE}-mutserve"

    bcftools reheader -s <(echo "${OLD_NAME} ${NEW_NAME}") ${VCF_FILE} \
    | bcftools +setGT -- -t a -n c:'1' \
    | bcftools annotate --rename-chrs ${chrom_map} -x INFO,^FORMAT/GT,^FORMAT/AF \
    | bcftools norm -f ${REF} --multiallelics -any \
    | bcftools norm -f ${REF} --atomize --atom-overlaps . \
    | bcftools filter -e 'FILTER != "PASS" && FMT/AF < 0.95' --mode x \
    | bcftools view -f PASS -t '^chrM:3107' -Oz -o ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz

    bcftools index --tbi ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz
    echo "${NEW_NAME},${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz" >> ${MERGE_MANIFEST}

done

## mitoscope pacbio
awk -F',' '{print $1,$2}' ${MITOSCOPE_PB_MANIFEST} | while read -r SAMPLE VCF_FILE
do  
    tech='pacbio'
    echo "mitoscope pacbio"
    echo $SAMPLE
    mkdir -p ${OUTDIR}/${tech}

    OLD_NAME=$(bcftools query -l "$VCF_FILE")
    NEW_NAME="${SAMPLE}-mitoscope"

    bcftools reheader -s <(echo "${OLD_NAME} ${NEW_NAME}") ${VCF_FILE} \
    | bcftools +setGT -- -t a -n c:'1' \
    | bcftools annotate --rename-chrs ${chrom_map} -x INFO,^FORMAT/GT,^FORMAT/HPL -c 'FORMAT/AF:=FORMAT/HPL' \
    | bcftools view -Oz -o ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz

    bcftools index --tbi ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz
    echo "${NEW_NAME},${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz" >> ${MERGE_MANIFEST}

done

## mitoscope ont
awk -F',' '{print $1,$2}' ${MITOSCOPE_ONT_MANIFEST} | while read -r SAMPLE VCF_FILE
do 
    tech='ont'
    echo "mitoscope ont"
    echo $SAMPLE
    mkdir -p ${OUTDIR}/${tech}

    OLD_NAME=$(bcftools query -l "$VCF_FILE")
    NEW_NAME="${SAMPLE}-mitoscope"

    bcftools reheader -s <(echo "${OLD_NAME} ${NEW_NAME}") ${VCF_FILE} \
    | bcftools +setGT -- -t a -n c:'1' \
    | bcftools annotate --rename-chrs ${chrom_map} -x INFO,^FORMAT/GT,^FORMAT/HPL -c 'FORMAT/AF:=FORMAT/HPL'  \
    | bcftools view -Oz -o ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz

    bcftools index --tbi ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz
    echo "${NEW_NAME},${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz" >> ${MERGE_MANIFEST}

done

## himito pacbio
awk -F',' '{print $1,$2}' ${HIMITO_PB_MANIFEST} | while read -r SAMPLE VCF_FILE
do
    tech='pacbio'
    echo "himito pacbio"
    echo $SAMPLE
    mkdir -p ${OUTDIR}/${tech}
    
    OLD_NAME=$(bcftools query -l "$VCF_FILE")
    NEW_NAME="${SAMPLE}-himito"

    bcftools reheader -s <(echo "${OLD_NAME} ${NEW_NAME}") ${VCF_FILE} \
    | bcftools +setGT -- -t a -n c:'1' \
    | bcftools annotate --rename-chrs ${chrom_map} -x INFO,^FORMAT/GT,^FORMAT/HF -c 'FORMAT/AF:=FORMAT/HF' \
    | bcftools view -f PASS -Oz -o ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz

    bcftools index --tbi ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz
    echo "${NEW_NAME},${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz" >> ${MERGE_MANIFEST}
done

## himito ont
awk -F',' '{print $1,$2}' ${HIMITO_ONT_MANIFEST} | while read -r SAMPLE VCF_FILE
do
    tech='ont'
    echo "himito ont"
    echo $SAMPLE
    mkdir -p ${OUTDIR}/${tech}
    
    OLD_NAME=$(bcftools query -l "$VCF_FILE")
    NEW_NAME="${SAMPLE}-himito"

    bcftools reheader -s <(echo "${OLD_NAME} ${NEW_NAME}") ${VCF_FILE} \
    | bcftools +setGT -- -t a -n c:'1' \
    | bcftools annotate --rename-chrs ${chrom_map} -x INFO,^FORMAT/GT,^FORMAT/HF -c 'FORMAT/AF:=FORMAT/HF' \
    | bcftools view -f PASS -Oz -o ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz

    bcftools index --tbi ${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz
    echo "${NEW_NAME},${OUTDIR}/${tech}/${NEW_NAME}.vcf.gz" >> ${MERGE_MANIFEST}
done

