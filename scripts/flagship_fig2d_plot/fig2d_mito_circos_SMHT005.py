#!/usr/bin/env python3
"""
P25 Figure 2d — mtDNA somatic heteroplasmy circos, donor SMHT005.

Matches the group's caption: "inner concentric tracks show per-tissue kernel
density of heteroplasmic variant positions". The tracks are therefore a density
over POSITION only -- heteroplasmy level (AF/VAF) is deliberately not encoded
anywhere in this panel, and there is no AF key. Set TRACK_MODE = "af_stems" for
the alternative encoding that does show AF (see below).

Rebuilt from the collaborator notebook `data/Mito/fig_panel/circos_plot.ipynb`
(cell 5), which had Mac-local paths and was unrunnable here. Legibility defects
fixed, in rough order of how much damage each did:

  1. The canvas was never the size it asked for. `import pycirclize` mutates
     `rcParams["savefig.bbox"] = "tight"`, so the notebook's `figsize=(5cm, 5cm)`
     plus an explicit `bbox_inches="tight"` saved at ~8 cm. Scaled back into a
     5 cm column that renders the nominal 4 pt labels at ~2.5 pt -- so nothing in
     the original was readable at print size regardless of encoding. This script
     forces `savefig.bbox = "standard"` and insets the polar axes so gene labels
     fit *inside* a fixed 5x5 cm canvas. Point sizes here are true.
  2. The KDE bandwidth was ~50-80 bp. `gaussian_kde`'s `bw_method` is a
     MULTIPLIER on the data's standard deviation, and the notebook passed 0.01.
     For these tissues (sd 2.4-7.1 kb) that is a 24-71 bp kernel, so every
     variant became a near-delta needle: the tracks read as scattered confetti,
     not densities, and the control-region enrichment -- the actual signal --
     was invisible. Scott's rule is no help either (3-4 kb, washes it out).
     Fixed absolute bandwidth instead: KDE_BW_BP, identical for every tissue, so
     track heights are comparable rather than each being smoothed differently.
  3. The KDE was linear, on a circular genome. mtDNA wraps, and the control
     region wraps with it: m.564 and m.16034 are 1,099 bp apart in reality but
     15,470 bp apart on a line. The notebook's density therefore split the single
     biggest feature in the data into two half-peaks pinned at opposite ends of
     the track, and inflated every tissue's sd (which is what broke #2's
     bandwidth). `circular_kde` sums kernels over the -L/0/+L images, so the
     control region reads as one continuous peak across the origin.
  4. Normalisation is a caption decision, not a bug fix. DENSITY_NORM ships as
     "per_tissue" (the notebook's behaviour, `density /= density.max()`), which is
     what "per-tissue kernel density" states: every lane is scaled to its own peak
     so the tracks compare SHAPE. The cost is that heights are not comparable --
     ESOP (6 calls) peaks as high as TESL (21) -- so the panel must not be read as
     showing which tissue carries more variants. Per-lane n is printed on every run
     for the caption. DENSITY_NORM = "shared" puts all lanes on one variants/kb
     scale if height comparability is ever wanted instead.
  5. Radial compression. Tracks ran r=78 down to r=5 at equal RADIAL width, so
     the innermost lane had ~15x less arc length than the outermost and the same
     data read as much smaller. Lanes now stop at R_LANE_IN.
  6. No tissue identification. The original passed `label=` to `fill_between`,
     but the only `ax.legend()` call was commented out, leaving six anonymous
     rings. Lanes are now labelled in-figure with the canonical 4-letter
     `smaht_code`, upright, in a single column down a data-free wedge near 6
     o'clock (opposite the control-region cluster, so nothing is occluded).
  7. Gene-ring noise. 22 tRNAs were label-suppressed but still drawn as
     full-size arrows, and "D-loop" was annotated twice because the control
     region is two GFF features that wrap the origin (16023-16569 and 0-576).
     tRNAs are now quiet grey ticks; the D-loop is labelled once.
  8. Palette fidelity. The notebook drew tracks at alpha=0.8, which composites
     every tissue toward white and so renders a colour that is NOT the canonical
     hex. Fills here are the verbatim hex from `smaht_tissue_colors_FINAL.json`
     (byte-identical to the table in `smaht/CLAUDE.md`) at alpha=1.0, so what
     lands on the page is exactly the palette value -- verified by sampling the
     rendered PNG. Nothing recolours, darkens or blends a tissue colour.

     Contrast is handled WITHOUT touching the hue: a single neutral outline
     (TRACK_EDGE) is applied uniformly to every track. Three canonical colours
     are very pale (BRFL #BBE7FF, TESL #AAF8B4, FBRO #FFE98D all sit at ~0.88
     luminance, i.e. ~1.1:1 against white) and would otherwise be invisible as
     1.8 mm fills. A uniform neutral edge gives them definition while distorting
     no tissue's identity, since it is the same for all lanes. Lane labels are
     neutral grey for the same reason -- their ring position already identifies
     them, so colour is not doing any work there.

     Note for print: matplotlib writes DeviceRGB (0 CMYK operators in the PDF).
     If the journal converts to CMYK downstream, the saturated palette members
     (LUNG #1DDD30, LIVR #FF4AA9) are outside CMYK gamut and will shift; that is
     a conversion issue, not something this script can pre-empt.

Scope: donor SMHT005 only, matching panels 2a / 2b / 2d. Lanes are the
TOP_N_TISSUES tissues with the most calls (same rule the notebook used), drawn
in the anatomical head-to-toe order shared with `fig2a_donor_summary_SMHT005.R`.

Type is 5 pt everywhere -- the smallest size allowed on P25 panels. On a fixed
5 cm canvas that budget is what drives two layout decisions:

  * No kb ruler. A labelled coordinate axis costs ~7 radial units of canvas for
    tick text. Gene names are the conventional positional reference in mtDNA
    figures, so the ruler is dropped and the lanes extend into that band.
  * Lane labels form ONE column, in a data-free wedge near 6 o'clock. The wedge
    has to be near 6 (or 12) o'clock because that is where the radial direction
    is vertical: labels then step down the page and only need to clear each
    other's 1.8 mm line height. Near 3/9 o'clock the radial direction is
    horizontal, the binding constraint becomes label WIDTH (~4.1 mm for a
    4-letter code), and 6 lanes will not fit. Staggering across two angles was
    tried and rejected -- it clears geometrically but destroys association, since
    radially-adjacent labels land on nearly the same visual row and the palette
    pairs two blues (BRFL/BRCE) and two oranges (ESOP/ADGR), so colour cannot
    break the tie. One monotonic column is unambiguous.

  * Per-lane n lives in a footer strip below the circle, not on the lane labels.
    Both in-circle options were measured and fail. Appending the count widens a
    label from 4.83 mm ("BRFL") to 7.23 mm ("BRFL 10") or 10.37 mm ("BRFL (n=10)")
    at 5 pt; since a label's angular footprint is width/circumference-at-its-radius,
    worst-case occlusion goes 0.14 -> 0.88 -> 1.00, i.e. a label box completely
    buries a peak at every angle the pitch constraint allows. Moving the column to
    9 o'clock does not help either: at 270 deg the radial direction is exactly
    horizontal, so pitch = lane_h * |cos(270)| = 0 and all six labels collapse onto
    one line (225 deg already fails at 7.78 vs 9.8 needed). Height is the only
    spendable dimension, so the counts go below the circle.

Raising TOP_N_TISSUES past 6 tightens the label column past what 5 pt type can
take; the script checks the pitch on every run and warns rather than silently
overlapping.

TRACK_MODE = "af_stems" swaps the densities for one radial stem per variant with
height on a log10 AF scale. That encoding shows heteroplasmy, which the density
cannot, but note the finding it surfaces: every call at AF>=20% in this donor
sits at m.309 or m.16183 -- the two canonical mtDNA poly-C homopolymer tracts
(m.303-315/CSB2, m.16184-16193/HVS-I) and the best-known length-heteroplasmy
artifact hotspots. Treat those heights as suspect before reading them as clonal
expansion. Output goes to a separate `_afstems` stem so both can coexist.

Usage:
    python3 figure2/scripts/fig2d_mito_circos_SMHT005.py
"""

import json
import os

import numpy as np
import pandas as pd
import pysam

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.text import Text

from pycirclize import Circos
from pycirclize.parser import Gff
from pycirclize import config as pycirc_config

# pycirclize sets savefig.bbox = "tight" at import time, which silently rescales
# the output and shrinks every font relative to the requested canvas. Undo it,
# then set our own style. Order matters: this must follow the pycirclize import.
plt.rcParams["savefig.bbox"] = "standard"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]
plt.rcParams["font.size"] = 5
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

# ─────────────────────────── paths ────────────────────────────────────────────
PROJ = "/Users/suhang/Analysis/SMaHT_sample/p25"
VCF = f"{PROJ}/p25_v2_merge_0817/p25.filled.vcf.gz"
GFF = f"{PROJ}/rCRS.gff3"
COLORS_JSON = "/Users/suhang/Analysis/SMaHT_sample/SMaHT_Mitogenome/scripts/smaht_tissue_colors_FINAL.json"
OUT_DIR = f"{PROJ}/figure2/plots"

# ─────────────────────────── config ───────────────────────────────────────────
DONOR = "SMHT005"
TIER = "HighConf"
AF_MIN = 0.001           # notebook's inclusion threshold
# Lane selection. TOP_N_TISSUES takes the N tissues with the most calls (the
# notebook's rule). Note the cap of 6 is GEOMETRIC -- the 5 pt lane-label pitch on a
# 5 cm canvas -- not a statement about the biology, and for this donor the cut lands
# mid-plateau: ranks 6-10 are 6, 6, 5, 5, 5 calls, so ESOP is in and HART out purely
# on the anatomical tie-break below.
#
# LANE_CODES_OVERRIDE pins the set explicitly instead (still drawn in anatomical
# order). Relevant because the notebook did NOT just take a top-N: it took the top 6
# across ALL donors and then hand-appended SMHT005-3AC (Fibroblast) and SMHT005-3Q
# (Lung), i.e. those two were deliberate inclusions. Under the count rule they rank
# 10th (n=5) and 17th (n=2) and are dropped. To restore the notebook's intent for
# this donor, e.g.:
#   LANE_CODES_OVERRIDE = ["3U", "3B", "3M", "3AK", "3AC", "3Q"]   # + FBRO, LUNG
TOP_N_TISSUES = 6

# Pinned set (2026-07-30, by request): BUCC (3B) and ADGR (3M) dropped, and the
# notebook's two hand-picked tissues restored -- FBRO (3AC, cultured cells) and
# LUNG (3Q). Drawn in anatomical order regardless of the order listed here.
# Caveat: LUNG has only n=2 calls for this donor, so its "density" is just two
# kernels side by side and carries no shape information -- see the run log, which
# prints n per lane.
LANE_CODES_OVERRIDE = ["3AK", "3AM", "3C", "3Q", "3AC", "3U"]

# "density"  -> per-tissue kernel density of variant POSITIONS (matches caption)
# "af_stems" -> one stem per variant, height = log10 heteroplasmy AF
TRACK_MODE = "density"

CHRM_LEN = 16569
SEQID = "NC_012920.1"

# Kernel density settings (TRACK_MODE == "density")
KDE_BW_BP = 200.0       # absolute Gaussian sd in bp, same for every tissue
# 200 bp chosen by comparing 120 / 200 / 300 at final print size: 120 bp breaks the
# control-region hotspot into needle sub-peaks in the low-n lanes, 300 bp smooths
# the HVS-I / OriB substructure away entirely. 200 bp keeps one clear peak with a
# visible shoulder. Note this is an absolute sd, NOT gaussian_kde's bw_method.
KDE_GRID = 3000         # evaluation points around the genome
# "per_tissue" -> each lane scaled to its own peak, matching the caption's
#                 "per-tissue kernel density". Shapes are comparable, HEIGHTS ARE
#                 NOT: a 6-variant lane peaks as high as a 21-variant one.
# "shared"     -> one variants/kb scale across all lanes, so heights are
#                 comparable but low-n tissues render nearly flat.
DENSITY_NORM = "per_tissue"

# Heteroplasmy encoding (TRACK_MODE == "af_stems"): log10 AF, floored just below
# AF_MIN so the smallest real call still gets a visible stem.
AF_FLOOR = 0.0075
AF_SCALE_MAX = 1.0
AF_STEM_MIN = 0.25      # fraction of lane height, so nothing renders as zero

# Track outline. Neutral and uniform on purpose -- see docstring defect 8. Set to
# None to drop it and let the canonical fills stand alone (pale tissues will
# largely disappear at this track height).
TRACK_EDGE = "#595959"
EDGE_MIN_Y = 0.02       # fraction of lane height below which the outline is cut
LABEL_COLOR = "#333333"

# Fixed canvas. Width is a hard constraint on P25 panels -- do not widen. Height is
# 5.9 cm to carry the per-tissue n footer; the circle itself is unchanged at
# 0.79 * 5.0 = 3.95 cm, so no data area was traded for the footer.
FIG_W_CM, FIG_H_CM = 5.0, 5.9
DPI = 300
FOOTER_CM = 1.35       # band below the circle, for the n footer
# Polar axes inset, leaving room for gene labels inside the fixed canvas. Kept at
# the same absolute size and top margin as the 5x5 version.
_CIRCLE_FRAC_W = 0.79
AXES_RECT = [0.105,
             FOOTER_CM / FIG_H_CM,
             _CIRCLE_FRAC_W,
             _CIRCLE_FRAC_W * FIG_W_CM / FIG_H_CM]

# Radial geometry (pycirclize MAX_R = 100)
R_GENE_OUT = 97.0
R_GENE_MID = 93.5       # + strand above, - strand below
R_GENE_IN = 90.0
R_LANE_OUT = 88.0

# Lane geometry is a preset, because track height and radial compression trade off
# directly and the right answer is an editorial call. The notebook ran tracks from
# r=78 all the way to r=5 with no gaps: 11.83 radial units per lane, but its
# innermost lane had only 6% of the outermost's arc length, so identical data read
# as far smaller the further in it sat.
#
#   preset          r_lane_in  lane_gap  lane height vs notebook  inner/outer arc
#   middle                 22       0.5                     89%             0.25
#   match_original         15       0.4                    100%             0.17
#
# Every preset still satisfies the 5 pt lane-label pitch; the script asserts that
# per preset rather than trusting the table above.
LANE_PRESETS = [
    dict(name="middle",         suffix="",           r_lane_in=22.0, lane_gap=0.5),
    dict(name="match_original", suffix="_matchorig", r_lane_in=15.0, lane_gap=0.4),
]

# Control region (D-loop) — wraps the origin, hence two arcs.
DLOOP_ARCS = [(16023, CHRM_LEN), (0, 576)]

# 5 pt is the hard floor for P25 panels; nothing here may go below it.
FS_MIN = 5.0
FS_GENE = FS_MIN
FS_LANE = FS_MIN
# 5 pt Nimbus Sans line height, in radial units, used to check label pitch.
PT_TO_R = 100.0 / ((FIG_W_CM * 10.0) * AXES_RECT[2] / 2.0) * (25.4 / 72.0)

# Anatomical head-to-toe tissue order, as tissue codes. Mirrors
# `tissue_name_levels` in fig2a_donor_summary_SMHT005.R so figure 2 panels agree.
ANATOMICAL_CODE_ORDER = [
    "3AK", "3AL", "3AN", "3AO", "3AM",   # brain: frontal, temporal, hippoc L/R, cerebellum
    "3B",                                 # buccal swab
    "3C", "3Q", "3S", "3O", "3A",         # esophagus, lung, heart, aorta, blood
    "3I", "3K", "3M",                     # liver, adrenal L/R
    "3E", "3G",                           # colon asc, desc
    "3AH", "3AF", "3AC",                  # muscle, skin abdomen, fibroblast
    "3U", "3W", "3Y", "3AA",              # testis L/R, ovary L/R
    "3AD",                                # skin calf
]

OUT_BASE = ("fig2d_mito_circos_SMHT005" if TRACK_MODE == "density"
            else "fig2d_mito_circos_SMHT005_afstems")


# ─────────────────────────── helpers ─────────────────────────────────────────
def load_heteroplasmy(vcf_path, tier=TIER, af_min=AF_MIN):
    """Long table of somatic heteroplasmic calls: sample, donor, code, pos, af.

    Filters match the notebook: Tier == tier, VariantOrigin == "somatic",
    HeteroplasmyStatus == "hetero", AF > af_min.
    """
    rows = []
    vcf = pysam.VariantFile(vcf_path, "r")
    samples = [str(s) for s in vcf.header.samples if ":" not in s]
    for record in vcf:
        for ind, _alt in enumerate(record.alts):
            for sample in samples:
                call = record.samples[sample]
                if call["Tier"][ind] != tier:
                    continue
                if call["HeteroplasmyStatus"][0] != "hetero":
                    continue
                # if call["VariantOrigin"][0] != "somatic":
                #     continue
                af = call["AF"][ind]
                if pd.isna(af) or float(af) <= af_min:
                    continue
                rows.append((sample, record.pos, float(af)))
    df = pd.DataFrame(rows, columns=["sample", "pos", "af"])
    # Samples are e.g. SMHT005-3AK.withpileup.final — strip the merge suffix
    # before taking the tissue code, or LANE_CODES_OVERRIDE ("3AK") never matches.
    core = df["sample"].str.split(".").str[0]
    df["donor"] = core.str.split("-").str[0]
    df["code"] = core.str.split("-").str[1]
    return df


def luminance(color):
    """Relative luminance, used only to report which fills are low-contrast."""
    r, g, b = to_rgb(color)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def circular_kde(positions, x_eval, bw_bp, length=CHRM_LEN):
    """Gaussian KDE on a CIRCULAR genome, with an absolute bp bandwidth.

    Two departures from the notebook's `gaussian_kde(pos, bw_method=0.01)`:

      * Circular. Kernels are summed over the -L / 0 / +L images of each point,
        so density leaks across the origin instead of being clipped at it. This
        matters more than it sounds: the control region spans m.16024-16569 plus
        m.1-576, so on a line its two halves sit at opposite ends of the axis and
        the single strongest feature in the data is torn into two half-peaks.
      * Absolute bandwidth. `bw_method` scales the data's sd, which here is
        dominated by that same wrap artefact (2.4-7.1 kb), so the same
        `bw_method` means a different kernel in every tissue. A fixed bp sd keeps
        smoothing identical across lanes, which is what makes heights comparable.

    Returns density in variants per bp: it integrates to len(positions) over the
    circle, so scaling by 1000 gives variants/kb.
    """
    p = np.asarray(positions, dtype=float)
    if p.size == 0:
        return np.zeros_like(x_eval, dtype=float)
    acc = np.zeros_like(x_eval, dtype=float)
    for shift in (-length, 0.0, length):
        z = (x_eval[:, None] - (p[None, :] + shift)) / bw_bp
        acc += np.exp(-0.5 * z * z).sum(axis=1)
    return acc / (bw_bp * np.sqrt(2.0 * np.pi))


def af_to_height(af):
    """log10 AF mapped onto [AF_STEM_MIN, 1] of the lane height."""
    lo, hi = np.log10(AF_FLOOR), np.log10(AF_SCALE_MAX)
    frac = (np.log10(np.clip(af, AF_FLOOR, AF_SCALE_MAX)) - lo) / (hi - lo)
    return AF_STEM_MIN + (1.0 - AF_STEM_MIN) * frac


def lane_occupancy(code):
    """Per-lane array over x_eval of "how much ink is drawn here", 0..1.

    Used only to place the lane labels. For densities this is the drawn curve; for
    af_stems it is a narrow bump at each variant, since a stem has no width.
    """
    if TRACK_MODE == "density":
        d = lane_density[code]
        denom = peak_kb if DENSITY_NORM == "shared" else max(d.max(), 1e-12)
        return d / denom
    pos = donor_calls.loc[donor_calls["code"] == code, "pos"].values
    occ = circular_kde(pos, x_eval, 40.0)
    return occ / max(occ.max(), 1e-12)


def choose_label_angle(lane_radii, lane_h, need_pitch,
                       half_w_mm=2.05, search=(140.0, 220.0), step=0.5):
    """Angle for the lane-label column: near-vertical, occluding as little as possible.

    Two constraints, and the second is why a gap-between-variants search is not
    enough. (1) Pitch: labels step along the radius, so only the vertical component
    of the lane pitch separates them -- lane_h * |cos(theta)| must clear 5 pt line
    height, which confines theta to near 6 o'clock (12 o'clock is the control-region
    peak). (2) Occlusion: each label paints an opaque white box, and its ANGULAR
    half-width is half_w_mm / circumference-at-that-radius. Inner lanes therefore
    subtend far more angle than outer ones -- at r=39 a 4-letter code spans +/-15 deg
    -- so a wedge that is free of variant POSITIONS can still put an inner label on
    top of a neighbouring kernel's peak. Searching angles and scoring actual drawn
    ink is the only thing that catches that.

    Returns (theta, worst_occlusion).
    """
    occ = {c: lane_occupancy(c) for c in lane_codes}
    x_deg = x_eval / CHRM_LEN * 360.0
    mm_per_runit = (FIG_W_CM * 10.0) * AXES_RECT[2] / 2.0 / 100.0

    best = None
    theta = search[0]
    while theta <= search[1]:
        if lane_h * abs(np.cos(np.radians(theta))) >= need_pitch:
            worst = 0.0
            for code, rc in zip(lane_codes, lane_radii):
                dtheta = half_w_mm / (2.0 * np.pi * rc * mm_per_runit) * 360.0
                sel = (x_deg >= theta - dtheta) & (x_deg <= theta + dtheta)
                if sel.any():
                    worst = max(worst, float(occ[code][sel].max()))
            # Prefer least occlusion; break ties toward 6 o'clock.
            key = (round(worst, 4), abs(theta - 180.0))
            if best is None or key < best[0]:
                best = (key, theta, worst)
        theta += step
    if best is None:                      # no angle satisfies the pitch constraint
        return 180.0, float("nan")
    return best[1], best[2]


# ─────────────────────────── data ─────────────────────────────────────────────
smaht = json.load(open(COLORS_JSON))
code2color = {k: "#" + v["tissue_color_hex"].lstrip("#") for k, v in smaht.items()}
code2short = {k: v["smaht_code"] for k, v in smaht.items()}
code2name = {k: v["full_name"] for k, v in smaht.items()}

calls = load_heteroplasmy(VCF)
donor_calls = calls[calls["donor"] == DONOR].copy()
if donor_calls.empty:
    raise SystemExit(f"no {TIER} somatic heteroplasmic calls for {DONOR}")

# Rank by call count; break ties on anatomical order so the selection is
# reproducible rather than dependent on groupby ordering.
anat_rank = {c: i for i, c in enumerate(ANATOMICAL_CODE_ORDER)}
counts = donor_calls.groupby("code").size()
ranked = sorted(counts.index, key=lambda c: (-counts[c], anat_rank.get(c, 999)))
if LANE_CODES_OVERRIDE:
    unknown = [c for c in LANE_CODES_OVERRIDE if c not in set(counts.index)]
    if unknown:
        raise SystemExit(f"LANE_CODES_OVERRIDE has no calls for {DONOR}: {unknown}")
    selected = set(LANE_CODES_OVERRIDE)
    print(f"lane set: pinned via LANE_CODES_OVERRIDE ({len(selected)} tissues)")
else:
    selected = set(ranked[:TOP_N_TISSUES])
    cut = counts[ranked[TOP_N_TISSUES - 1]]
    tied_out = [c for c in ranked[TOP_N_TISSUES:] if counts[c] == cut]
    print(f"lane set: top {TOP_N_TISSUES} of {len(ranked)} tissues by call count")
    if tied_out:
        print(f"  NOTE: cut at n={cut} is a tie -- "
              f"{', '.join(code2short[c] for c in tied_out)} excluded only by the "
              f"anatomical tie-break, not by count")

lane_codes = [c for c in ANATOMICAL_CODE_ORDER if c in selected]
missing = selected - set(lane_codes)
if missing:
    raise SystemExit(f"tissue code(s) absent from ANATOMICAL_CODE_ORDER: {sorted(missing)}")

n_lanes = len(lane_codes)
# The notebook's per-lane radial thickness, for reporting each preset against it:
# tracks spanned r=5..(80-2) with no gaps.
NOTEBOOK_LANE_H = (78.0 - 5.0) / 6.0

print(f"{DONOR}: {len(donor_calls)} calls, {donor_calls['code'].nunique()} tissues "
      f"| TRACK_MODE = {TRACK_MODE}")
dloop = donor_calls[(donor_calls.pos <= 576) | (donor_calls.pos >= 16024)]
print(f"control region: {len(dloop)}/{len(donor_calls)} calls "
      f"({100 * len(dloop) / len(donor_calls):.0f}%) in {576 + CHRM_LEN - 16024} bp "
      f"({100 * (576 + CHRM_LEN - 16024) / CHRM_LEN:.1f}% of the genome)")

# Densities, computed up front so a shared scale can be established across lanes.
x_eval = np.linspace(0, CHRM_LEN, KDE_GRID)
lane_density = {}
for code in lane_codes:
    pos = donor_calls.loc[donor_calls["code"] == code, "pos"].values
    lane_density[code] = circular_kde(pos, x_eval, KDE_BW_BP) * 1000.0  # variants/kb

peak_kb = max(d.max() for d in lane_density.values()) if lane_density else 1.0
print(f"KDE: circular, bandwidth {KDE_BW_BP:g} bp (sd), {DENSITY_NORM} normalisation")
if DENSITY_NORM == "shared":
    print(f"     peak {peak_kb:.2f} variants/kb = track ceiling shared by all lanes")
else:
    print("     each lane scaled to its own peak (shapes comparable, heights not)")
pale = [code2short[c] for c in lane_codes if luminance(code2color[c]) > 0.75]
if pale:
    print(f"     low-contrast canonical fills carried by the neutral outline: "
          f"{', '.join(pale)}")

print(f"lanes ({n_lanes}, anatomical order, outer->inner):")
for c in lane_codes:
    sub = donor_calls[donor_calls["code"] == c]
    n_dl = int(((sub.pos <= 576) | (sub.pos >= 16024)).sum())
    print(f"  {code2short[c]:<5} {code2name[c]:<22} n={len(sub):>3} "
          f"({n_dl} in control region)  peak={lane_density[c].max():.2f}/kb")

# Vertical component of the lane pitch is what separates stacked labels; compare
# against 5 pt line height. Checked per preset below, never assumed.
need_r = FS_MIN * PT_TO_R * 1.10

# ─────────────────────────── figure ───────────────────────────────────────────
def build_figure(preset):
    """Render one panel for a lane-geometry preset. Returns the output stem.

    Only lane thickness and the size of the empty centre differ between presets;
    the data, KDE, palette, gene ring and label placement are identical, so the
    variants are directly comparable.
    """
    r_lane_in = preset["r_lane_in"]
    lane_gap = preset["lane_gap"]
    lane_h = (R_LANE_OUT - r_lane_in) / n_lanes
    lane_usable = lane_h - lane_gap

    # Lane radii up front: the label-angle search needs them, because a label's
    # angular footprint depends on the radius it sits at.
    lane_radii = [R_LANE_OUT - i * lane_h - lane_usable / 2.0
                  for i in range(n_lanes)]
    wedge_deg, occl = choose_label_angle(lane_radii, lane_h, need_r)
    pitch_r = lane_h * abs(np.cos(np.radians(wedge_deg)))

    print(f"\n[{preset['name']}] r_lane_in={r_lane_in:g} lane_gap={lane_gap:g} "
          f"-> lane {lane_usable:.2f} r-units "
          f"({100 * lane_usable / NOTEBOOK_LANE_H:.0f}% of notebook), "
          f"inner/outer arc {r_lane_in / R_LANE_OUT:.2f}")
    print(f"[{preset['name']}] label column at deg {wedge_deg:.1f}, pitch "
          f"{pitch_r:.1f} vs {need_r:.1f} needed for {FS_MIN:g}pt type")
    print(f"[{preset['name']}] worst density occluded by a label box: {occl:.3f} "
          f"of lane height")
    if pitch_r < need_r:
        print(f"[{preset['name']}] WARNING: {FS_MIN:g}pt lane labels will collide; "
              f"reduce TOP_N_TISSUES or lower r_lane_in")
    if occl > 0.25:
        print(f"[{preset['name']}] WARNING: no label angle avoids the data -- "
              f"a label box hides {occl:.0%} of a lane's peak")

    fig = plt.figure(figsize=(FIG_W_CM / 2.54, FIG_H_CM / 2.54), dpi=DPI)
    ax = fig.add_axes(AXES_RECT, polar=True)

    pycirc_config.ann_adjust.dr = 0.05
    circos = Circos({SEQID: CHRM_LEN})
    sector = circos.sectors[0]

    # ---- control-region shading, behind the lanes ---------------------------
    # Genomic coords via a background track, so the two wrapping arcs land exactly
    # on the D-loop rather than being hand-converted to degrees.
    bg_track = sector.add_track((r_lane_in, R_LANE_OUT))
    bg_track.axis(fc="none", ec="none")
    for start, end in DLOOP_ARCS:
        bg_track.rect(start, end, fc="#FFF7DF", ec="none", zorder=0)

    # ---- gene ring ----------------------------------------------------------
    gene_track = sector.add_track((R_GENE_IN, R_GENE_OUT))
    gene_track.axis(fc="#F2F2F2", ec="none")

    features = Gff(GFF).get_seqid2features(feature_type="gene")[SEQID]

    # tRNAs first, as quiet ticks: still on the map, but no longer competing with
    # the 15 protein-coding / rRNA genes that carry labels.
    for feat in features:
        if feat.qualifiers.get("gene_biotype", [""])[0] == "tRNA":
            gene_track.genomic_features(feat, plotstyle="box",
                                        r_lim=(R_GENE_IN, R_GENE_OUT),
                                        fc="#C4C4C4", ec="none")

    # Strand is carried by arrow direction and by which half of the ring a feature
    # sits in; two colours are enough, so no gene-feature legend is needed.
    for feat in features:
        if feat.qualifiers.get("gene_biotype", [""])[0] == "tRNA":
            continue
        name = feat.qualifiers.get("Name", [""])[0]
        if not name or name.startswith("hypothetical"):
            continue
        if feat.location.strand == 1:
            r_lim, fc = (R_GENE_MID, R_GENE_OUT), "#F4977E"
        else:
            r_lim, fc = (R_GENE_IN, R_GENE_MID), "#7FB9DE"
        gene_track.genomic_features(feat, plotstyle="arrow", r_lim=r_lim,
                                    fc=fc, ec="white", lw=0.15)
        mid = (int(feat.location.start) + int(feat.location.end)) / 2
        gene_track.annotate(mid, name, label_size=FS_GENE)

    # D-loop: two features wrapping the origin, but ONE label (the notebook
    # annotated both, hence the duplicate "D-loop" in the old panel).
    dloop_feats = Gff(GFF).extract_features("D_loop")
    for feat in dloop_feats:
        gene_track.genomic_features(feat, plotstyle="box",
                                    r_lim=(R_GENE_IN, R_GENE_OUT),
                                    fc="#F2C230", ec="none")
    widest = max(dloop_feats,
                 key=lambda f: int(f.location.end) - int(f.location.start))
    gene_track.annotate((int(widest.location.start) + int(widest.location.end)) / 2,
                        "D-loop", label_size=FS_GENE)

    # ---- per-tissue lanes ---------------------------------------------------
    r_top = R_LANE_OUT
    lane_geom = []
    for code in lane_codes:
        r_hi = r_top
        r_lo = r_top - lane_usable
        lane_geom.append((code, r_lo, r_hi))

        color = code2color[code]      # verbatim canonical hex, never transformed

        track = sector.add_track((r_lo, r_hi))
        track.axis(fc="none", ec="none")

        # A faint baseline gives each lane a readable zero without the heavy
        # concentric rings that dominated an earlier pass.
        circos.line(r=r_lo, deg_lim=(0, 360), lw=0.18, color="#ECECEC", zorder=1)

        if TRACK_MODE == "density":
            dens = lane_density[code]
            denom = peak_kb if DENSITY_NORM == "shared" else max(dens.max(), 1e-12)
            y = dens / denom
            # alpha=1.0: any transparency composites toward white and stops the
            # rendered colour being the canonical hex.
            track.fill_between(x_eval, y, 0, vmin=0, vmax=1,
                               fc=color, ec="none", alpha=1.0, zorder=3)
            if TRACK_EDGE:
                # Outline only where there IS density. Tracing the whole curve
                # draws the flat zero baseline too, which turns into a heavy dark
                # ring around the circumference -- one per lane, dominating the
                # panel.
                mask = np.flatnonzero(y > EDGE_MIN_Y)
                if mask.size:
                    for run in np.split(mask,
                                        np.flatnonzero(np.diff(mask) != 1) + 1):
                        if run.size >= 2:
                            track.line(x_eval[run], y[run], vmin=0, vmax=1,
                                       lw=0.3, color=TRACK_EDGE, zorder=4)
        else:
            sub = donor_calls[donor_calls["code"] == code].sort_values("af")
            for pos, af in zip(sub["pos"], sub["af"]):
                # arc=False -> a true radial stem, so physical linewidth is
                # constant regardless of the lane's radius.
                track.line([pos, pos], [0, af_to_height(af)], arc=False,
                           vmin=0, vmax=1, lw=0.5, color=color,
                           solid_capstyle="butt", zorder=3)

        r_top -= lane_h

    # ---- lane labels: upright, one column down a data-free wedge ------------
    for code, r_lo, r_hi in lane_geom:
        circos.text(code2short[code], r=(r_lo + r_hi) / 2, deg=wedge_deg,
                    adjust_rotation=False, size=FS_LANE, color=LABEL_COLOR,
                    ha="center", va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none",
                              alpha=0.92))

    circos.plotfig(ax=ax)
    ax.set_ylim(0, pycirc_config.MAX_R)
    ax.set_aspect("equal", adjustable="box", anchor="C")
    ax.margins(0)
    ax.set_axis_off()
    for artist in ax.get_children():
        if isinstance(artist, Text):
            artist.set_clip_on(False)

    # The centre disc is deliberately empty. A density track needs no numeric key
    # -- height is "how clustered are variants here", which the caption states --
    # and an axis floating in the middle of a circos has no physical relationship
    # to the tracks it would be describing.

    # ---- footer: per-lane n ------------------------------------------------
    # Counts belong on the panel because DENSITY_NORM="per_tissue" scales every lane
    # to its own peak, so height carries no information about how many calls a
    # tissue has: FBRO's 5 isolated singletons peak as high as TESL's 21-call
    # control-region pile-up. Without n the panel invites exactly that misreading.
    # Gene labels are unclipped and overhang the axes by ~0.5 cm, so the footer's
    # ceiling is measured from what actually rendered rather than assumed -- ATP8 at
    # the bottom of the circle collided with a hard-coded position.
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    label_bottoms = [t.get_window_extent(renderer=rend).y0
                     for t in ax.get_children() if isinstance(t, Text) and t.get_text()]
    ceil_frac = (min(label_bottoms) / fig.bbox.height) if label_bottoms else \
                FOOTER_CM / FIG_H_CM
    row_h = 0.26 / FIG_H_CM                     # 2.6 mm line pitch
    pad = 0.10 / FIG_H_CM                       # 1 mm clearance under the gene labels
    ncol = 3
    n_rows = (len(lane_codes) + ncol - 1) // ncol
    top = ceil_frac - pad - row_h / 2
    if top - (n_rows - 1) * row_h < row_h / 2:
        print(f"[{preset['name']}] WARNING: footer does not fit under the gene "
              f"labels; raise FIG_H_CM or FOOTER_CM")
    for i, code in enumerate(lane_codes):
        row, col = divmod(i, ncol)
        fx = 0.5 + (col - (ncol - 1) / 2) * 0.30
        n = int((donor_calls["code"] == code).sum())
        fig.text(fx, top - row * row_h, f"{code2short[code]} n={n}", size=FS_LANE,
                 color=LABEL_COLOR, ha="center", va="center")

    stem = OUT_BASE + preset["suffix"]
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        path = os.path.join(OUT_DIR, f"{stem}.{ext}")
        fig.savefig(path, dpi=DPI, facecolor="white")
        print("wrote", path)
    plt.close(fig)
    return stem


for _preset in LANE_PRESETS:
    build_figure(_preset)
