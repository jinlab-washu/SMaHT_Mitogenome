#!/usr/bin/env python3
"""
Transform a single-tissue VCF from p25_v2:
  - Moves the INFO_TO_FORMAT fields from INFO to FORMAT (see FORMAT_SPECS)
  - Collapses per-caller sample columns into one per-tissue sample column
    (GT: diploid 0/1 if any caller called the alt, 0/0 if ref, ./. if missing;
     AF: mean of non-missing caller AF)
  - Writes uncompressed VCF to stdout (pipe to bgzip for indexing)

Usage:
    python3 transform_vcf.py input.vcf.gz SAMPLE_NAME | bgzip > output.vcf.gz
    tabix -p vcf output.vcf.gz
"""

import gzip
import sys

# INFO fields to move into FORMAT. Number/Type/Description copied from p25_v2
# (*.withpileup.final.vcf.gz) INFO headers so merge keeps one value per sample.
FORMAT_SPECS = [
    ("Tier", "A", "String", "Confidence tier: HighConf,LowConf,LikelyArtifact"),
    ("TissueSupport", "A", "String", "Variant support across tissues"),
    ("AvailableTech", "A", "String", "List of available technologies used for sequencing."),
    ("AvailableTissues", "A", "String", "List of available tissues sequenced from donor."),
    ("HeteroplasmyStatus", "A", "String", "Heteroplasmy status (homoplasmic, heteroplasmic)"),
    ("VariantOrigin", "A", "String", "Variant origin (germline, somatic)"),
    ("TechSupport_Caller", ".", "String", "Variant support across technologies (illumina, pacbio, ont) based on caller support."),
    ("TechSupport_Pileup", ".", "String", "Variant support across technologies (illumina, pacbio, ont) based on mpileup data support."),
    ("TechSupport_Combined", ".", "String", "Variant support across technologies (illumina, pacbio, ont) based on combined caller and mpileup data support."),
    ("CallerSupport", ".", "String", "Variant support across callers (mutect2, mutserve, mitoscope, himito)"),
    ("OriginalTier", ".", "String", "The original confidence tier from caller before pileup incorporation."),
    ("Illumina_AF_caller", "1", "Float", "Averaged heteroplasmy frequency across Illumina callers."),
    ("PacBio_AF_caller", "1", "Float", "Averaged heteroplasmy frequency across PacBio callers."),
    ("ONT_AF_caller", "1", "Float", "Averaged heteroplasmy frequency across ONT callers."),
    ("Illumina_AF_pileup", "1", "Float", "Variant allele fraction based on Illumina pileup."),
    ("PacBio_AF_pileup", "1", "Float", "Variant allele fraction based on PacBio pileup."),
    ("ONT_AF_pileup", "1", "Float", "Variant allele fraction based on ONT pileup."),
]

INFO_TO_FORMAT = [name for name, *_ in FORMAT_SPECS]
NEW_FORMAT_HEADERS = [
    f'##FORMAT=<ID={name},Number={number},Type={typ},Description="{desc}">'
    for name, number, typ, desc in FORMAT_SPECS
]


def parse_info(info_str):
    d = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            d[k] = v
        else:
            # Flag-only INFO (e.g. TechSupport_Pileup with no techs): keep
            # as missing so FORMAT stays '.' rather than the Python True.
            d[field] = '.'
    return d


def parse_sample(fmt_str, sample_str):
    keys = fmt_str.split(':')
    vals = sample_str.split(':')
    # pad with '.' if sample has fewer fields than format keys
    vals += ['.'] * (len(keys) - len(vals))
    return dict(zip(keys, vals))


def rebuild_info(info_str, skip_keys):
    """Return INFO string with the specified keys removed."""
    parts = []
    for field in info_str.split(';'):
        key = field.split('=')[0] if '=' in field else field
        if key not in skip_keys:
            parts.append(field)
    return ';'.join(parts) if parts else '.'


def transform(input_path, sample_name, out=sys.stdout):
    opener = gzip.open if input_path.endswith('.gz') else open
    with opener(input_path, 'rt') as fh:
        format_headers_inserted = False
        for line in fh:
            line = line.rstrip('\n')

            # ── meta-information lines ──────────────────────────────────────
            if line.startswith('##'):
                # Drop INFO definitions for fields we are moving to FORMAT
                if any(f'##INFO=<ID={f},' in line for f in INFO_TO_FORMAT):
                    continue
                # Standard diploid GT header (inputs sometimes use Number=A / haploid)
                if line.startswith('##FORMAT=<ID=GT,'):
                    print(
                        '##FORMAT=<ID=GT,Number=1,Type=String,'
                        'Description="Genotype">',
                        file=out,
                    )
                    continue
                print(line, file=out)

            # ── #CHROM header ───────────────────────────────────────────────
            elif line.startswith('#CHROM'):
                # Insert new FORMAT headers just before the column header line
                for h in NEW_FORMAT_HEADERS:
                    print(h, file=out)
                cols = line.split('\t')
                # Replace all sample columns with the single tissue name
                print('\t'.join(cols[:9] + [sample_name]), file=out)

            # ── data lines ──────────────────────────────────────────────────
            else:
                cols = line.split('\t')
                chrom, pos, id_, ref, alt, qual, filt, info_str, fmt_str = cols[:9]
                caller_cols = cols[9:]

                info = parse_info(info_str)
                callers = [parse_sample(fmt_str, c) for c in caller_cols]

                # ── collapse GT across callers → haploid (mtDNA) ──────────
                # mtDNA is haploid: emit 1 (alt), 0 (ref), . (missing) so
                # fill-tags gives AN = one allele/sample and AF = carrier fraction.
                # (Diploid 0/1 would double AN and cap AF near 0.5.)
                def _has_alt(g):
                    if g in ('.', './.', '.|.', '', None):
                        return False
                    return any(a not in ('0', '.') for a in g.replace('|', '/').split('/'))

                def _is_missing(g):
                    return g in ('.', './.', '.|.', '', None)

                gt_vals = [c.get('GT', '.') for c in callers]
                if any(_has_alt(g) for g in gt_vals):
                    gt = '1'
                elif all(_is_missing(g) for g in gt_vals):
                    gt = '.'
                else:
                    gt = '0'

                # ── collapse per-caller AF (average of non-missing) ─────────
                af_vals = [float(c['AF']) for c in callers if c.get('AF', '.') != '.']
                af = f"{sum(af_vals) / len(af_vals):.4f}" if af_vals else '.'

                # ── INFO fields to promote ──────────────────────────────────
                # Build FORMAT/sample from INFO_TO_FORMAT so the two stay in
                # sync: any tag added to that list is promoted automatically.
                promoted_vals = [info.get(f, '.') for f in INFO_TO_FORMAT]

                new_fmt    = ':'.join(['GT', 'AF'] + INFO_TO_FORMAT)
                new_sample = ':'.join([gt, af] + promoted_vals)

                new_info = rebuild_info(info_str, set(INFO_TO_FORMAT))

                print('\t'.join([chrom, pos, id_, ref, alt, qual, filt,
                                 new_info, new_fmt, new_sample]), file=out)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: {sys.argv[0]} input.vcf.gz SAMPLE_NAME', file=sys.stderr)
        sys.exit(1)
    transform(sys.argv[1], sys.argv[2])
