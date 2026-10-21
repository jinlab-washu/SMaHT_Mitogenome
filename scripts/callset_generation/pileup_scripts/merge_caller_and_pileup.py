#!/usr/bin/env python3

# ─────────────────────────────────────────────────────────────────────────────
# Pileup support is incorporated to update existing confidence tiers determined
# based on caller support. Variant annotation fields including Tier are updated 
# and additional fields are added to reflect incorporated pileup results 
# (TechSupport_Pileup, TechSupport_Caller, TechSupport_Combined). Heteroplasmy
# levels from pileup data are reported separately alongside caller derived values
# individually for each technology. Tier status is updated from LikelyArtifact 
# or LowConf to HighConf if additional tech support is detected using pileup data
# resulting in variant support from at least two different seq technologies
# across caller and pileup derived support.
# ─────────────────────────────────────────────────────────────────────────────

import argparse
from pysam import VariantFile
import sys
import math

def get_args():
    parser = argparse.ArgumentParser(description="Incorporate pileup support to update existing variant annotation fields such as tier, tech support, and AF.")
    parser.add_argument("-s", "--sample", help="Sample name used to subset minipileup VCFs.", required=True)
    parser.add_argument("-c", "--caller_input", help="File path of integrated caller VCF (subset to one donor).", required=True)
    parser.add_argument("-p", "--minipileup_pacbio", help="File path of minipileup Pacbio VCF.", required=True)
    parser.add_argument("-n", "--minipileup_ont", help="File path of minipileup ONT VCF.", required=True)
    parser.add_argument("-i", "--minipileup_illumina", help="File path of minipileup Illumina VCF.", required=True)
    parser.add_argument("-o", "--output", help="File path for output VCF.", required=True)
    return parser.parse_args()
args = get_args()

donor_tissue = args.sample
fr_caller = VariantFile(args.caller_input, 'r')
fr_mini_pb = VariantFile(args.minipileup_pacbio, 'r')
fr_mini_ont = VariantFile(args.minipileup_ont, 'r')
fr_mini_illumina = VariantFile(args.minipileup_illumina, 'r')


fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "TechSupport_Caller"),
        ("Number", "."),
        ("Type", "String"),
        ("Description", "Variant support across technologies (illumina, pacbio, ont) based on caller support."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "TechSupport_Pileup"),
        ("Number", "."),
        ("Type", "String"),
        ("Description", "Variant support across technologies (illumina, pacbio, ont) based on mpileup data support."),
    ],
)
fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "TechSupport_Combined"),
        ("Number", "."),
        ("Type", "String"),
        ("Description", "Variant support across technologies (illumina, pacbio, ont) based on combined caller and mpileup data support."),
    ],
)
fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "CallerSupport"),
        ("Number", "."),
        ("Type", "String"),
        ("Description", "Variant support across callers (mutect2, mutserve, mitoscope, himito)"),
    ],
)
fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "OriginalTier"),
        ("Number", "."),
        ("Type", "String"),
        ("Description", "The original confidence tier from caller before pileup incorporation."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "PacBio_AF_caller"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Averaged heteroplasmy frequency across PacBio callers."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "ONT_AF_caller"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Averaged heteroplasmy frequency across ONT callers."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "Illumina_AF_caller"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Averaged heteroplasmy frequency across Illumina callers."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "PacBio_AF_pileup"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Variant allele fraction based on PacBio pileup."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "ONT_AF_pileup"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Variant allele fraction based on ONT pileup."),
    ],
)

fr_caller.header.add_meta(
    "INFO",
    items=[
        ("ID", "Illumina_AF_pileup"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Variant allele fraction based on Illumina pileup."),
    ],
)

fw = VariantFile(args.output, 'w', header=fr_caller.header)

class SampleVariant:
    """Store variant information for a single variant from pileup data (position, ref, alt, VAF)."""
    def __init__(self, id, vaf):
        self.id = id
        self.vaf = vaf

    
def populate_minipileup_pass_variants(fr_mini, mini_pass_variants):
    """
    Read in minipileup variants that pass filtration criteria (PassAltReadSupport sites with SBP > 0.01 or 0.3 < SBR < 0.7).
    Populations mini_pass_variants dictionary with sample name as key and list of SampleVariant objects as value.
    
    """
    for record in fr_mini.fetch():
        ref = record.ref
        alt = record.alts[0]
        pos = record.pos

        for sample_name in record.samples:
            record_sample = record.samples[sample_name]
            if record_sample['RS'] == 'PassAltReadSupport' and (record_sample['SBP'] >= 0.01 or (record_sample['SBR'] >= 0.3 and record_sample['SBR'] <= 0.7)):
                variant_id=f'{pos}:{ref}>{alt}'

                if sample_name in mini_pass_variants:
                    mini_pass_variants[sample_name].append(SampleVariant(variant_id, record_sample['AF']))
                else:
                    mini_pass_variants[sample_name] = [SampleVariant(variant_id, record_sample['AF'])]

def find_variant_in_minipileup(variant_id, sample_name, minipileup_dict):
    """
    Find variant in minipileup dict using sample and variant ID (ie 3243:A>G).
    Returns VAF of identified variant if present, otherwise None.
    """
    if sample_name in minipileup_dict:
        target_variant = next((v for v in minipileup_dict[sample_name] if v.id == variant_id), None)
        if target_variant:
            return target_variant.vaf
        else:
            return None
    else:
        return None

mini_pb_pass_variants = {}
mini_ont_pass_variants = {}
mini_illumina_pass_variants = {}

### population dict for each tech with minipileup variants that pass filtration criteria (RS, SBP, SBR)
populate_minipileup_pass_variants(fr_mini_pb, mini_pb_pass_variants)
populate_minipileup_pass_variants(fr_mini_ont, mini_ont_pass_variants)
populate_minipileup_pass_variants(fr_mini_illumina, mini_illumina_pass_variants)

#desired_order = ["Tier", "TechSupport", "PileupSupport", "ToolSupport", "TissueSupport", 'AvailableTech', 'AvailableTissues', 'Illumina_AF', 'PacBio_AF', 'ONT_AF', 'HeteroplasmyStatus', 'VariantOrigin', 'NewTier']

for record in fr_caller.fetch():
    pileupsupport = []

    ref = record.ref
    alt = record.alts[0]
    pos = record.pos
    variant_id=f'{pos}:{ref}>{alt}'

    tier = record.info['Tier'][0]
    techsupport = sorted(record.info['TechSupport'][0].split('+'))
    toolsupport = sorted(record.info['ToolSupport'][0].split('+'))
    tissuesupport = sorted(record.info['TissueSupport'][0].split('+'))

    pb_var_af = find_variant_in_minipileup(variant_id, f'{donor_tissue}-pacbio', mini_pb_pass_variants)
    if pb_var_af:
        pileupsupport.append('pacbio')

    ont_var_af = find_variant_in_minipileup(variant_id, f'{donor_tissue}-ont', mini_ont_pass_variants)
    if ont_var_af:
        pileupsupport.append('ont')

    illumina_var_af = find_variant_in_minipileup(variant_id, f'{donor_tissue}-illumina', mini_illumina_pass_variants)
    if illumina_var_af:
        pileupsupport.append('illumina')

    ## get list of pileup tech support + generate new tech support list based on combined caller + pileup support
    pileupsupport = sorted(pileupsupport)
    combtechsupport = sorted(set(techsupport + pileupsupport))
    
    ## add new INFO fields for tech support based on pileup, caller and combined
    record.info['TechSupport_Combined'] = '+'.join(combtechsupport)
    record.info['TechSupport_Pileup'] = '+'.join(pileupsupport)
    record.info['TechSupport_Caller'] = '+'.join(techsupport)
    record.info['CallerSupport'] = '+'.join(toolsupport)

    ## store old tier and update final tier based on whether additional tech support from pileup detected
    record.info['OriginalTier'] = tier
    if len(combtechsupport) > 1:
        record.info['Tier'] = 'HighConf'
        if set(combtechsupport) != set(techsupport):
            print(f'{variant_id},{tier},{techsupport},{combtechsupport},{pileupsupport},{toolsupport},{tissuesupport}')
    else:
        record.info['Tier'] = tier

    ## add in VAF info for new variant support from pileup (retain original VAF if already made by caller)?
    record.info['Illumina_AF_caller'] = record.info.get('Illumina_AF') 
    record.info['PacBio_AF_caller'] = record.info.get('PacBio_AF')
    record.info['ONT_AF_caller'] = record.info.get('ONT_AF')

    record.info['Illumina_AF_pileup'] = illumina_var_af
    record.info['PacBio_AF_pileup'] = pb_var_af
    record.info['ONT_AF_pileup'] = ont_var_af
    
    fw.write(record)


fr_caller.close()
fr_mini_pb.close()
fr_mini_ont.close()
fr_mini_illumina.close()
fw.close()