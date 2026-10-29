"""
Divergent gap figure using Hi-C compartment delta (INTS11 vs ctl).
Single bigWig, unstranded. Same layout as divergent_gap_figure.py:
  profile (mean ± SEM) on top, heatmap below, pairs sorted by gap size.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import pyBigWig
from collections import defaultdict

# ── Parameters ────────────────────────────────────────────────────────────────
MIN_GAP    = 1_000
MAX_GAP    = 10_000
FLANK_BP   = 3_000
GAP_BINS   = 100
FLANK_BINS = 60
TOTAL_BINS = FLANK_BINS + GAP_BINS + FLANK_BINS
TSS_B_BIN  = FLANK_BINS
TSS_A_BIN  = FLANK_BINS + GAP_BINS

BW_PATH = "data/compartments/delta_INT11vsctl.bw"

cmap = mcolors.LinearSegmentedColormap.from_list(
    "blue_white_orange",
    ["#4393C3", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#E08214"],
    N=512,
)

# ── Load BED and find divergent pairs ─────────────────────────────────────────
def load_bed(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip(): continue
            cols = line.strip().split('\t')
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name = cols[3]; strand = cols[5]
            tss = start if strand == '+' else end
            genes.append(dict(chrom=chrom, start=start, end=end,
                              name=name, strand=strand, tss=tss))
    return genes

def find_divergent_pairs(genes, min_gap, max_gap):
    by_chrom = defaultdict(list)
    for g in genes:
        by_chrom[g['chrom']].append(g)
    for v in by_chrom.values():
        v.sort(key=lambda g: g['tss'])
    pairs = []
    for cg in by_chrom.values():
        for i, ga in enumerate(cg):
            if ga['strand'] != '+': continue
            for j in range(i - 1, -1, -1):
                gb = cg[j]
                if gb['tss'] >= ga['tss']: continue
                if gb['strand'] == '-':
                    gap = ga['tss'] - gb['tss']
                    if min_gap <= gap <= max_gap:
                        pairs.append((gb, ga, gap))
                break
    return pairs

# ── Extract scaled profile for one pair ───────────────────────────────────────
def extract_profile(bw, chrom_lengths, gb, ga):
    chrom_bw = "chr" + ga['chrom']
    clen  = chrom_lengths.get(chrom_bw)
    out   = np.full(TOTAL_BINS, np.nan, dtype=np.float32)
    if clen is None:
        return out
    tss_b = gb['tss']
    tss_a = ga['tss']

    def fetch(s, e, n):
        s, e = max(0, s), min(clen, e)
        if s >= e: return None
        return bw.stats(chrom_bw, s, e, type='mean', nBins=n)

    def fill(offset, vals):
        if vals:
            for k, v in enumerate(vals):
                if v is not None:
                    out[offset + k] = float(v)

    fill(0,                     fetch(tss_b - FLANK_BP, tss_b, FLANK_BINS))
    fill(FLANK_BINS,            fetch(tss_b, tss_a, GAP_BINS))
    fill(FLANK_BINS + GAP_BINS, fetch(tss_a, tss_a + FLANK_BP, FLANK_BINS))
    return out

# ── Main ──────────────────────────────────────────────────────────────────────
genes = load_bed("genes_01_bydistance.bed")
pairs = find_divergent_pairs(genes, MIN_GAP, MAX_GAP)
pairs.sort(key=lambda p: p[2])   # sort by gap size ascending
n_pairs = len(pairs)
print(f"Divergent pairs (gap {MIN_GAP//1000}–{MAX_GAP//1000} kb): {n_pairs}")

bw     = pyBigWig.open(BW_PATH)
clens  = dict(bw.chroms())

mat = np.full((n_pairs, TOTAL_BINS), np.nan, dtype=np.float32)
for i, (gb, ga, gap) in enumerate(pairs):
    mat[i] = extract_profile(bw, clens, gb, ga)
bw.close()

prof = np.nanmean(mat, axis=0)
sem  = np.nanstd(mat, axis=0) / np.sqrt(n_pairs)
x    = np.arange(TOTAL_BINS)

# symmetric color limits at 99th percentile of abs values
vlim = np.nanpercentile(np.abs(mat), 99)

# ── x-axis ticks ─────────────────────────────────────────────────────────────
x_ticks  = [0, FLANK_BINS // 2, TSS_B_BIN,
             TSS_A_BIN, TSS_A_BIN + FLANK_BINS // 2, TOTAL_BINS]
x_labels = ["-3kb", "-1.5kb", "TSS_B\n(−)", "TSS_A\n(+)", "+1.5kb", "+3kb"]

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(9, 10))
gs  = gridspec.GridSpec(2, 1, height_ratios=[1.5, 5], hspace=0.08)

# ─ Profile ───────────────────────────────────────────────────────────────────
ax_prof = fig.add_subplot(gs[0])
ax_prof.fill_between(x, prof - sem, prof + sem, color="#4393C3", alpha=0.25)
ax_prof.plot(x, prof, color="#4393C3", linewidth=1.2)
ax_prof.axhline(0, color="grey", linewidth=0.6, linestyle="--")
ax_prof.axvline(TSS_B_BIN, color="black", linewidth=0.8, linestyle="--")
ax_prof.axvline(TSS_A_BIN, color="black", linewidth=0.8, linestyle="--")
ax_prof.axvspan(TSS_B_BIN, TSS_A_BIN, color="grey", alpha=0.08)
ax_prof.set_xlim(0, TOTAL_BINS)
ax_prof.set_xticks([])
ax_prof.set_ylabel("Mean Δ compartment\n(INTS11 vs ctl)", fontsize=8)
ax_prof.set_title(
    f"Hi-C compartment delta across divergent gene pairs\n"
    f"n = {n_pairs} pairs, gap {MIN_GAP//1000}–{MAX_GAP//1000} kb  |  "
    f"shaded = gap region  |  mean ± SEM",
    fontsize=9, fontweight="bold"
)

# ─ Heatmap ───────────────────────────────────────────────────────────────────
ax_heat = fig.add_subplot(gs[1])
im = ax_heat.imshow(
    mat, aspect="auto", cmap=cmap,
    vmin=-vlim, vmax=vlim,
    interpolation="nearest",
    extent=[0, TOTAL_BINS, n_pairs, 0],
)
ax_heat.axvline(TSS_B_BIN, color="black", linewidth=0.8, linestyle="--")
ax_heat.axvline(TSS_A_BIN, color="black", linewidth=0.8, linestyle="--")
ax_heat.set_xlim(0, TOTAL_BINS)
ax_heat.set_xticks(x_ticks)
ax_heat.set_xticklabels(x_labels, fontsize=7)
ax_heat.set_ylabel("Gene pairs\n(gap size ↑)", fontsize=8)
ax_heat.tick_params(labelsize=7)
ax_heat.set_title("Δ compartment score  (INTS11 vs ctl)", fontsize=9,
                  fontweight="bold", loc="left")

# gap size axis on right
ax2 = ax_heat.twinx()
gaps = [p[2] for p in pairs]
ax2.set_ylim(n_pairs, 0)
ax2.set_ylabel("Gap size (bp)", fontsize=7)
ax2.set_yticks(np.linspace(0, n_pairs, 5))
ax2.set_yticklabels(
    [f"{int(g):,}" for g in np.interp(np.linspace(0, n_pairs, 5),
                                       [0, n_pairs], [gaps[0], gaps[-1]])],
    fontsize=7
)

cbar = plt.colorbar(im, ax=ax_heat, shrink=0.6, pad=0.12,
                    label="Δ compartment score")
cbar.ax.tick_params(labelsize=7)

plt.savefig("divergent_gap_hic.svg", bbox_inches="tight")
plt.savefig("divergent_gap_hic.png", bbox_inches="tight", dpi=150)
print("Saved divergent_gap_hic.svg / .png")
