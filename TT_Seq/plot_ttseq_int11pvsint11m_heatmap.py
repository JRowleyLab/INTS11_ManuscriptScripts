"""
Heatmap + profile of TT-seq log2 fold-change (INTS11+ vs INTS11−) around TSS.
Same regions / order as the Hi-C delta heatmap:
genes in genes_01_bydistance.bed order, 20 kb upstream / 4 kb downstream.
Signal: data/ttseq_bigwigs/log2_int11pvsint11m.bw
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors

UPSTREAM   = 20000
DOWNSTREAM = 4000
BINSIZE    = 50
N_BINS     = (UPSTREAM + DOWNSTREAM) // BINSIZE   # 480
TSS_BIN    = UPSTREAM // BINSIZE                   # 400

cmap = mcolors.LinearSegmentedColormap.from_list(
    "blue_white_orange",
    ["#4393C3", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#E08214"],
    N=512,
)

# ── load matrix ───────────────────────────────────────────────────────────────
d   = np.load("ttseq_int11pvsint11m_stranded.npz", allow_pickle=True)
mat = d["matrix"]          # [8874 × 480]
n_genes = mat.shape[0]
print(f"Matrix: {mat.shape[0]:,} genes × {mat.shape[1]} bins")

# ── profile ───────────────────────────────────────────────────────────────────
prof = np.nanmean(mat, axis=0)
sem  = np.nanstd(mat, axis=0) / np.sqrt(n_genes)
x    = np.arange(N_BINS)

# fixed color limits
vmin, vmax = -0.4, 0.4
print(f"Color limits: {vmin} to {vmax}")

# ── x-axis ticks ─────────────────────────────────────────────────────────────
x_ticks  = [0, 100, 200, 300, TSS_BIN, TSS_BIN + 40, TSS_BIN + 80]
x_labels = ["-20kb", "-15kb", "-10kb", "-5kb", "TSS", "+2kb", "+4kb"]

# ── figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(8, 12))
gs  = gridspec.GridSpec(2, 1, height_ratios=[1.2, 5], hspace=0.06)

# profile
ax_p = fig.add_subplot(gs[0])
ax_p.fill_between(x, 0, prof, color="darkred", alpha=0.7)
ax_p.fill_between(x, prof - sem, prof + sem, color="darkred", alpha=0.25)
ax_p.plot(x, prof, color="darkred", linewidth=1.0)
ax_p.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax_p.axvline(TSS_BIN, color="black", linewidth=0.8)
ax_p.set_xlim(0, N_BINS)
ax_p.set_xticks([])
ax_p.set_ylabel("Mean TT-seq log2FC\n(INTS11+ vs INTS11-)", fontsize=8)
ax_p.set_title(
    "TT-seq log2 fold-change (INTS11+ vs INTS11-) around TSS\n"
    f"n = {n_genes:,} genes  |  sorted by upstream gene distance  |  mean ± SEM",
    fontsize=9, fontweight="bold"
)

# heatmap
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
                    label="TT-seq log2FC (INTS11+ vs INTS11-)")
cbar.ax.tick_params(labelsize=7)

plt.savefig("ttseq_delta_INT11pvsINT11m_heatmap.svg", bbox_inches="tight")
plt.savefig("ttseq_delta_INT11pvsINT11m_heatmap.png", bbox_inches="tight", dpi=150)
print("Saved ttseq_delta_INT11pvsINT11m_heatmap.svg / .png")
