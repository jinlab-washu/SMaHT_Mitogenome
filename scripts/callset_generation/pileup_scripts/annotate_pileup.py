#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
# This script parses a VCF from minipileup and adds additional FORMAT annotation
# fields such as coverage-based alt read support and strand bias checks. 
# ─────────────────────────────────────────────────────────────────────────────

import argparse
from pysam import VariantFile
import sys
import math
from scipy.stats import fisher_exact
from scipy.stats import poisson

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="File path of minipileup VCF.", required=True)
    parser.add_argument("-o", "--output", help="File path for output VCF.", required=True)
    parser.add_argument("-e", "--error_rate", help="Estimated sequencing error rate.", type=float, default=0.001, required=False)
    return parser.parse_args()
args = get_args()

error_rate = args.error_rate 
input_path= args.input
output_path= args.output

def get_read_cutoff(cov: int, error_rate: float = error_rate, target_p: float = 1e-2):
    """
    Compute read count cutoffs using Poisson approximation to the binomial distribution.
    Modified from SMaHT SNV working group pipeline (https://github.com/smaht-dac/calling-pipelines/blob/main/dockerfiles/calling_utils/scripts/tier_filter_variants_SR_PB_ONT.py)

    Parameters
    ----------
    cov : int
        Coverage for variant (ref + alt read depth)
    error_rate : float, optional
        Sequencing error rate (default = 0.001)
    target_p : float, optional
        Target probability for random error (default = 0.01)
    """

    def find_cutoff(cov: int, error_rate: float, target_p: float):
        """
        Find minimal count where Poisson tail < target_p.
        Modified from SMaHT SNV working group pipeline (https://github.com/smaht-dac/calling-pipelines/blob/main/dockerfiles/calling_utils/scripts/tier_filter_variants_SR_PB_ONT.py)
        """

        lmbda = cov * error_rate
        ## assume r will be at least lmda so can skip checking 1...lmda
        start = max(1, int(lmbda))

        for r in range(start, cov + 1):
            pt = poisson.sf(r - 1, lmbda)
            if pt < target_p:
                return r
        return cov

    cutoff = find_cutoff(cov, error_rate, target_p)

    return cutoff

def fisher_strand_bias(REF_ADF, REF_ADR, ALT_ADF, ALT_ADR, AF):
    """
    Compute Fisher's exact test p-value for strand bias using counts from minipileup.
    Modified from SMaHT SNV working group pipeline (https://github.com/smaht-dac/calling-pipelines/blob/main/dockerfiles/calling_utils/scripts/tier_filter_variants_SR_PB_ONT.py)
    """
    
    # Fisher's exact test
    table = [[REF_ADF, REF_ADR],
             [ALT_ADF, ALT_ADR]]

    ## for variants with high VAF, skip fisher test
    if AF > 0.75:
        p_value = 1.0
    else:
        try:
            _, p_value = fisher_exact(table, alternative="two-sided")
        except Exception as e:
            print(f"Fisher test failed for table {table}")
            p_value = 1.0

    ratio = ALT_ADF / (ALT_ADF + ALT_ADR)
    
    return (ratio,p_value)


fr = VariantFile(input_path, 'r')
fr.header.add_meta(
    "FORMAT",
    items=[
        ("ID", "SBP"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Strand bias p-value from Fisher Exact test."),
    ],
)
fr.header.add_meta(
    "FORMAT",
    items=[
        ("ID", "SBR"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Strand bias ratio based on alternate allele forward-strand fraction, calculated as ALT_ADF/(ALT_ADF+ALT_ADR); 0.5 indicates equal forward and reverse strand support"),
    ],
)
fr.header.add_meta(
    "FORMAT",
    items=[
        ("ID", "RS"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Alternate read count meets the minimum coverage-dependent Poisson error threshold [PassAltReadSupport, LowAltReadSupport, NoAltReadSupport]."),
    ],
)

fw = VariantFile(output_path, 'w', header=fr.header)


for record in fr.fetch():
    print(record.pos)
    for sample_name in record.samples:
        record_sample = record.samples[sample_name]

        gt = record_sample['GT']
        
        if gt == (None,None):
            record_sample['RS'] = None
            record_sample['SBP'] = None
            record_sample['SBR'] = None
        else:
            total_reads = int(record_sample['DPF'] + record_sample['DPR'])
            alt_reads = record_sample['ADF'][1] + record_sample['ADR'][1]
            ref_reads = record_sample['ADF'][0] + record_sample['ADR'][0]
            af = record_sample['AF']

            t = get_read_cutoff(total_reads)

            ## variant should have >= 2 read support on both F and R
            if alt_reads > t and record_sample['ADF'][1] >= 2 and record_sample['ADR'][1] >= 2:
                #print(dict(record_sample))
                rs = 'PassAltReadSupport'
            elif alt_reads == 0:
                rs = 'NoAltReadSupport'
            else:
                rs = 'LowAltReadSupport'

            ## calculate strand bias pval using fisher text + overall ratio of alt alelle forward vs total alt 
            if rs != 'NoAltReadSupport':
                sb_ratio, sb_pval = fisher_strand_bias(record_sample['DPF'] - record_sample['ADF'][1], record_sample['DPR'] - record_sample['ADR'][1], record_sample['ADF'][1], record_sample['ADR'][1], af)
            else:
                sb_ratio = 1.0
                sb_pval = 1.0

            record_sample['RS'] = rs
            record_sample['SBP'] = sb_pval
            record_sample['SBR'] = sb_ratio
        
    fw.write(record)


fr.close()
fw.close()