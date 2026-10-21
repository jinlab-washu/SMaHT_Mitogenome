#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
# Set up table with list of VCFs to merge per sample and the output VCF filepaths
# ─────────────────────────────────────────────────────────────────────────────

import pandas as pd
import sys

WORKDIR = sys.argv[1]
MANIFEST=f"{WORKDIR}/merge_manifests/p25.txt"
OUTFILE=f"{WORKDIR}/merge_manifests/p25.merge.txt"
OUTDIR=f"{WORKDIR}/output/merged"

df = pd.read_csv(MANIFEST, header=0, names=['Sample', 'FilePath'])
print(len(df))

df[['Donor', 'Tissue_Code', 'Tech', 'Tool']] = df['Sample'].str.strip().str.split('-', expand=True)

df_wide = df.groupby(['Donor', 'Tissue_Code'])['FilePath'].apply(list).reset_index()
df_wide['Donor_Tissue'] = df_wide['Donor'] + "-" + df_wide['Tissue_Code']

df_wide['FilePath'] = df_wide['FilePath'].apply(lambda x: ' '.join([item for item in x]))
df_wide['outpath'] = OUTDIR + "/" + df_wide['Donor_Tissue'] + ".vcf.gz"


print(df_wide)

df_wide.to_csv(OUTFILE, index=False)

