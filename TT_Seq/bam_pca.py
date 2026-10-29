import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

labels  = ["PAR_IAA_minus", "PAR_IAA_plus", "INTS1_minus", "INTS11_minus", "INTS11_plus"]
colors  = ["#2166AC", "#92C5DE", "#4DAC26", "#D6604D", "#F4A582"]
markers = ["o", "o", "s", "^", "^"]
PSEUDO  = 1.0   # added before log2; same units as CPM

def load_cpm(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    mat  = data["matrix"]
    m    = mat[~np.isnan(mat).any(axis=1)]
    totals = m.sum(axis=0)
    cpm  = m / (totals / 1e6)
    return cpm, totals

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle(
    "PCA — BAM CPM normalization (per million reads within BED regions)\n"
    "log₂(CPM + 1) transformed",
    fontsize=11, fontweight="bold"
)

datasets = [
    ("Promoters (TSS ±500 bp)",   "promoter_bam_summary.npz",  axes[0]),
    ("Gene Body + 1 kb upstream", "genebody_bam_summary.npz",  axes[1]),
]

for title, npz, ax in datasets:
    cpm, totals = load_cpm(npz)
    log_cpm = np.log2(cpm + PSEUDO)   # (regions, samples)

    # PCA: rows = samples, cols = regions
    pca    = PCA(n_components=2)
    coords = pca.fit_transform(log_cpm.T)   # (5, 2)
    pct    = pca.explained_variance_ratio_ * 100

    print(f"\n{'='*55}")
    print(f"{title}")
    print(f"  PC1 range: {coords[:,0].min():.4f}  to  {coords[:,0].max():.4f}  "
          f"  (span {coords[:,0].ptp():.4f})")
    print(f"  PC2 range: {coords[:,1].min():.4f}  to  {coords[:,1].max():.4f}  "
          f"  (span {coords[:,1].ptp():.4f})")
    print(f"  Variance: PC1={pct[0]:.1f}%  PC2={pct[1]:.1f}%")
    for i, l in enumerate(labels):
        print(f"    {l:<22} PC1={coords[i,0]:+.4f}  PC2={coords[i,1]:+.4f}")

    for i, (lbl, col, mk) in enumerate(zip(labels, colors, markers)):
        ax.scatter(coords[i, 0], coords[i, 1],
                   color=col, marker=mk, s=160, zorder=3,
                   edgecolors="black", linewidths=0.8, label=lbl)
        ax.annotate(lbl.replace("_", "\n"),
                    (coords[i, 0], coords[i, 1]),
                    textcoords="offset points", xytext=(8, 4), fontsize=7.5)

    # reference lines through origin
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--", zorder=0)
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--", zorder=0)

    # ---- axis range check ----
    # pad 20% around the data; do NOT let the zero-crossing force a huge axis
    x_pad = coords[:,0].ptp() * 0.20
    y_pad = coords[:,1].ptp() * 0.20
    ax.set_xlim(coords[:,0].min() - x_pad, coords[:,0].max() + x_pad)
    ax.set_ylim(coords[:,1].min() - y_pad, coords[:,1].max() + y_pad)

    ax.set_xlabel(f"PC1  ({pct[0]:.1f}% variance)", fontsize=10)
    ax.set_ylabel(f"PC2  ({pct[1]:.1f}% variance)", fontsize=10)
    ax.set_title(title, fontsize=10)

    # annotate library sizes in corner
    lib_str = "\n".join(f"{l}: {t/1e6:.1f}M" for l, t in zip(labels, totals))
    ax.text(0.02, 0.02, lib_str, transform=ax.transAxes, fontsize=6,
            va="bottom", color="#555555",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

plt.tight_layout()
plt.savefig("bam_pca.svg", bbox_inches="tight")
print("\nSaved bam_pca.svg")
