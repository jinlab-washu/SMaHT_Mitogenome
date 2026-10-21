#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Generate merged VCFs per donor-tissue sample using bcftools and simplified headers
# ─────────────────────────────────────────────────────────────────────────────


module load modules modules-init modules-gs
module load bcftools/1.21 htslib/1.21

WORKDIR=$1
MERGE_MANIFEST="${WORKDIR}/merge_manifests/p25.merge.txt"
TMPDIR="${WORKDIR}/output/tmp"
HEADER="${WORKDIR}/scripts/resources/header_for_tissue_merge.txt"
REF="${WORKDIR}/scripts/resources/MT.fasta"

mkdir -p ${TMPDIR}
mkdir -p ${WORKDIR}/output/merged

while IFS=',' read -r DONOR TISSUE_CODE LIST_OF_VCFS DONOR_TISSUE OUT_VCF
do
    echo $DONOR_TISSUE

    bcftools merge -m none --force-single ${LIST_OF_VCFS} \
    | bcftools norm -f ${REF} --multiallelics -any \
    | bcftools norm -f ${REF} --atomize --atom-overlaps . \
    | bcftools annotate -x 'INFO,FILTER,QUAL' \
    | bcftools +fill-tags -- -t AC \
    | bcftools view -t '^chrM:3107' -Oz -o ${OUT_VCF}.tmp

    cat ${HEADER} > ${TMPDIR}/header.tmp
    echo -e "$(zcat ${OUT_VCF}.tmp | grep '^#CHROM')" >> ${TMPDIR}/header.tmp
    bcftools reheader -h ${TMPDIR}/header.tmp ${OUT_VCF}.tmp -o ${OUT_VCF}
    
    bcftools index --tbi ${OUT_VCF} 

    ## rm temp files
    rm ${OUT_VCF}.tmp
    rm ${TMPDIR}/header.tmp

done < <(tail -n +2 ${MERGE_MANIFEST})

rmdir ${TMPDIR}
