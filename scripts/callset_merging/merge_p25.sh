#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# merge_p25.sh
#
# Pipeline: move per-sample INFO annotations (Tier, TissueSupport,
# AvailableTech/Tissues, HeteroplasmyStatus, VariantOrigin, CallerSupport,
# OriginalTier, TechSupport_{Caller,Pileup,Combined}, and per-tech
# Illumina/PacBio/ONT AF from callers and pileup) from INFO → FORMAT,
# collapse per-caller columns into one per-tissue sample, then merge all
# samples into a single multi-sample VCF.
#
# Requirements: python3, bgzip, tabix, bcftools
#
# Usage:
#   bash merge_p25.sh [INPUT_DIR] [OUTPUT_DIR] [THREADS]
#
# Optional env overrides:
#   REF_FASTA      default: /Users/suhang/Analysis/rCRS.fasta
#   OUTPUT_PREFIX  default: p25
#
# Defaults:
#   INPUT_DIR  = ../p25_v2
#   OUTPUT_DIR = ../p25_merged
#   THREADS    = 4
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${1:-${SCRIPT_DIR}/../p25_v2}"
OUTPUT_DIR="${2:-${SCRIPT_DIR}/../p25_merged}"
THREADS="${3:-4}"
REF_FASTA="${REF_FASTA:-/Users/suhang/Analysis/rCRS.fasta}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-p25}"
TRANSFORM="${SCRIPT_DIR}/transform_vcf.py"
# bcftools writes scratch files under TMPDIR; keep them off /tmp (sandbox/full-disk).
export TMPDIR="${TMPDIR:-${OUTPUT_DIR}/.tmp}"

[[ -f "${REF_FASTA}" ]] || { echo "missing REF_FASTA: ${REF_FASTA}" >&2; exit 1; }
[[ -f "${TRANSFORM}" ]] || { echo "missing transform script: ${TRANSFORM}" >&2; exit 1; }

mkdir -p "${OUTPUT_DIR}/per_sample" "${TMPDIR}"

echo "[1/4] Transforming VCFs (INFO → FORMAT, collapsing callers)..."
echo "  INPUT_DIR=${INPUT_DIR}"
echo "  OUTPUT_DIR=${OUTPUT_DIR}"
echo "  REF_FASTA=${REF_FASTA}"
echo "  THREADS=${THREADS}"

TRANSFORMED_LIST="${OUTPUT_DIR}/transformed.list"
: > "${TRANSFORMED_LIST}"

for vcf in "${INPUT_DIR}"/*.vcf.gz; do
    [ -f "${vcf}" ] || continue
    # SMHT001-3A.withpileup.final.vcf.gz → SMHT001-3A
    base="$(basename "${vcf}" .vcf.gz)"
    base="${base%.withpileup.final}"
    norm_out="${OUTPUT_DIR}/per_sample/${base}.norm.vcf.gz"

    if [ -f "${norm_out}" ] && [ -f "${norm_out}.tbi" ]; then
        echo "  skipping ${base} (already exists)"
    else
        echo "  ${base}"
        tmp_out="${OUTPUT_DIR}/per_sample/${base}.tmp.vcf.gz"
        # Split multiallelics + left-align indels to canonical biallelic records.
        python3 "${TRANSFORM}" "${vcf}" "${base}" \
            | bgzip -@ "${THREADS}" -c > "${tmp_out}"
        bcftools norm -f "${REF_FASTA}" -m -any "${tmp_out}" -Oz -o "${norm_out}"
        tabix -f -p vcf "${norm_out}"
        rm -f "${tmp_out}"
    fi

    echo "${norm_out}" >> "${TRANSFORMED_LIST}"
done

n_samples="$(wc -l < "${TRANSFORMED_LIST}" | tr -d ' ')"
echo "[2/4] Merging ${n_samples} samples with bcftools merge..."

MERGED_UNSORTED="${OUTPUT_DIR}/${OUTPUT_PREFIX}_merged_unsorted.vcf.gz"
# --merge none: keep every ALT as its own record (one full-sample row per allele)
# instead of packing alleles into a multiallelic diploid GT. mtDNA heteroplasmy can
# put 3+ alt alleles at one POS (e.g. pos 72 T>C/T>A/T>G); a diploid GT caps at 2 and
# would silently drop calls. Per-allele records have no such cap. Exact-duplicate
# records are removed later with `norm -d exact` (NOT -d both, which drops real alleles).
bcftools merge \
    --file-list "${TRANSFORMED_LIST}" \
    --merge none \
    --missing-to-ref \
    --output-type z \
    --threads "${THREADS}" \
    --output "${MERGED_UNSORTED}"

tabix -f -p vcf "${MERGED_UNSORTED}"

echo "[3/4] Sorting + split to biallelic + drop duplicate records..."

MERGED_SORTED="${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz"
MERGED_SORTED_TMP="${OUTPUT_DIR}/${OUTPUT_PREFIX}.sorted.tmp.vcf.gz"
bcftools sort \
    --output-type z \
    --output "${MERGED_SORTED_TMP}" \
    "${MERGED_UNSORTED}"
# -m -any: split multiallelic sites into biallelic records (one ALT per row)
# -d exact: drop only records identical in POS+REF+ALT (NOT -d both, which would
# collapse two distinct biallelic SNPs at the same POS and drop a real allele).
bcftools norm -f "${REF_FASTA}" -m -any -d exact -Oz -o "${MERGED_SORTED}" "${MERGED_SORTED_TMP}"
tabix -f -p vcf "${MERGED_SORTED}"
rm -f "${MERGED_UNSORTED}" "${MERGED_UNSORTED}.tbi" \
      "${MERGED_SORTED_TMP}" "${MERGED_SORTED_TMP}.tbi"

echo "[4/4] fill-tags..."

MERGED_FINAL="${OUTPUT_DIR}/${OUTPUT_PREFIX}.filled.vcf.gz"
bcftools +fill-tags "${MERGED_SORTED}" -O z -o "${MERGED_FINAL}" -- -t AF,AN,AC
tabix -f -p vcf "${MERGED_FINAL}"

echo ""
echo "Done! Output: ${MERGED_FINAL}"
echo ""
echo "Quick sanity check:"
bcftools stats "${MERGED_FINAL}" | grep -E "^SN"
