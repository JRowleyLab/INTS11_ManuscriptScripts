"""
Heatmap + profile of TT-seq log2 fold-change (INTS11+ vs INTS1−) around TSS.
Same regions / order as the Hi-C delta heatmap:
genes in genes_01_bydistance.bed order, 20 kb upstream / 4 kb downstream.
Signal: data/ttseq_bigwigs/log2_int11pvsint1m.bw

Input is a .npz matrix written by HiC/preprocessing/ttseq/stranded_computematrix.py
(one row per gene, one column per bin, TSS-anchored). Rows are plotted in the
order they appear in the .npz, which is the order of the BED file used to build it.

Running with no arguments reproduces the manuscript figure:
    python plot_ttseq_delta_heatmap.py

Running on other data:
    python plot_ttseq_delta_heatmap.py --npz my_matrix.npz --out-prefix my_heatmap \
        --vmin -1 --vmax 1 --title-label "Treatment vs Control"
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")   # draw to files only, no display needed (works on headless servers)
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors


def parse_args():
    """Read command-line options. Defaults reproduce the manuscript figure exactly."""
    p = argparse.ArgumentParser(description="TSS-anchored heatmap + mean profile from a stranded_computematrix.py .npz")
    p.add_argument("--npz", default="ttseq_int11pvsint1m_stranded.npz",
                   help="Input matrix from stranded_computematrix.py (default: %(default)s)")
    p.add_argument("--out-prefix", default="ttseq_delta_INT11vsctl_heatmap",
                   help="Output file prefix; writes <prefix>.svg and <prefix>.png (default: %(default)s)")
    p.add_argument("--vmin", type=float, default=-0.4, help="Lower colour limit (default: %(default)s)")
    p.add_argument("--vmax", type=float, default=0.4, help="Upper colour limit (default: %(default)s)")
    p.add_argument("--title-label", default="INTS11 vs INTS1",
                   help="Comparison name used in the title, y-axis and colourbar (default: %(default)s)")
    return p.parse_args()


def kb_label(bp):
    """Turn a signed bp offset from the TSS into a tick label like '-20kb', 'TSS' or '+2kb'."""
    if bp == 0:
        return "TSS"
    kb = bp / 1000
    # drop the trailing '.0' for whole kb values
    kb_str = f"{kb:g}"
    return f"+{kb_str}kb" if bp > 0 else f"{kb_str}kb"


args = parse_args()

cmap = mcolors.LinearSegmentedColormap.from_list(
    "blue_white_orange",
    ["#4393C3", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#E08214"],
    N=512,
)

# ── load matrix ───────────────────────────────────────────────────────────────
d   = np.load(args.npz, allow_pickle=True)
mat = d["matrix"]          # [n_genes × n_bins]; manuscript run: [8874 × 480]
n_genes = mat.shape[0]
print(f"Matrix: {mat.shape[0]:,} genes × {mat.shape[1]} bins")

# Window geometry is stored in the .npz by stranded_computematrix.py, so the
# figure adapts to whatever window/bin size the matrix was built with.
# Manuscript values: 20 kb upstream, 4 kb downstream, 50 bp bins → 480 bins, TSS at bin 400.
UPSTREAM   = int(d["upstream"])
DOWNSTREAM = int(d["downstream"])
BINSIZE    = int(d["bin_size"])
N_BINS     = (UPSTREAM + DOWNSTREAM) // BINSIZE
TSS_BIN    = UPSTREAM // BINSIZE

# ── profile ───────────────────────────────────────────────────────────────────
# Mean across genes for each bin, ignoring NaN (bins with no bigWig data).
prof = np.nanmean(mat, axis=0)
sem  = np.nanstd(mat, axis=0) / np.sqrt(n_genes)
x    = np.arange(N_BINS)

# fixed color limits
vmin, vmax = args.vmin, args.vmax
print(f"Color limits: {vmin} to {vmax}")

# ── x-axis ticks ─────────────────────────────────────────────────────────────
# Upstream: a tick every 5 kb back from the TSS. Downstream: a tick every 2 kb.
# For the manuscript window this gives -20,-15,-10,-5 kb, TSS, +2, +4 kb.
up_offsets   = list(range(-UPSTREAM, 0, 5000)) if UPSTREAM >= 5000 else [-UPSTREAM]
down_offsets = list(range(2000, DOWNSTREAM + 1, 2000)) if DOWNSTREAM >= 2000 else [DOWNSTREAM]
offsets  = up_offsets + [0] + down_offsets
x_ticks  = [TSS_BIN + off // BINSIZE for off in offsets]   # bp offset → bin index
x_labels = [kb_label(off) for off in offsets]

# ── figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(8, 12))
gs  = gridspec.GridSpec(2, 1, height_ratios=[1.2, 5], hspace=0.06)

# profile (top panel): mean ± SEM across all genes
ax_p = fig.add_subplot(gs[0])
ax_p.fill_between(x, 0, prof, color="darkred", alpha=0.7)
ax_p.fill_between(x, prof - sem, prof + sem, color="darkred", alpha=0.25)
ax_p.plot(x, prof, color="darkred", linewidth=1.0)
ax_p.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax_p.axvline(TSS_BIN, color="black", linewidth=0.8)
ax_p.set_xlim(0, N_BINS)
ax_p.set_xticks([])
ax_p.set_ylabel(f"Mean TT-seq log2FC\n({args.title_label})", fontsize=8)
ax_p.set_title(
    f"TT-seq log2 fold-change ({args.title_label}) around TSS\n"
    f"n = {n_genes:,} genes  |  sorted by upstream gene distance  |  mean ± SEM",
    fontsize=9, fontweight="bold"
)

# heatmap (bottom panel): one row per gene, in input order
ax_h = fig.add_subplot(gs[1])
im = ax_h.imshow(
    mat, aspect="auto", cmap=cmap,
    vmin=vmin, vmax=vmax,
    interpolation="nearest",
    extent=[0, N_BINS, n_genes, 0],
)
ax_h.axvline(TSS_BIN, color="black", linewidth=0.8)

ax_h.set_xticks(x_ticks)
ax_h.set_xticklabels(x_labels, fontsize=8)
ax_h.set_xlim(0, N_BINS)
ax_h.set_xlabel("Distance from TSS", fontsize=8)
ax_h.set_ylabel("Genes (closest upstream neighbor → top)", fontsize=8)
ax_h.tick_params(labelsize=7)

cbar = fig.colorbar(im, ax=ax_h, shrink=0.5, pad=0.02,
                    label=f"TT-seq log2FC ({args.title_label})")
cbar.ax.tick_params(labelsize=7)

plt.savefig(f"{args.out_prefix}.svg", bbox_inches="tight")
plt.savefig(f"{args.out_prefix}.png", bbox_inches="tight", dpi=150)
print(f"Saved {args.out_prefix}.svg / .png")
