import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

labels   = ["PAR_IAA_minus", "PAR_IAA_plus", "INTS1_minus", "INTS11_minus", "INTS11_plus"]
n        = len(labels)
PSEUDO   = 1e-2
MIN_EXPR = 1.0   # CPM

def load_cpm(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    mat  = data["matrix"]
    m    = mat[~np.isnan(mat).any(axis=1)]
    totals = m.sum(axis=0)
    return m / (totals / 1e6)

def build_median_abslfc(cpm):
    result = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = cpm[:, i], cpm[:, j]
            expressed = (a >= MIN_EXPR) | (b >= MIN_EXPR)
            lfc = np.log2((a[expressed] + PSEUDO) / (b[expressed] + PSEUDO))
            result[i, j] = np.median(np.abs(lfc))
    return result

def draw(ax, data_mat, vmax, title, n_genes):
    masked = np.where(np.eye(n, dtype=bool), np.nan, data_mat)
    im = ax.imshow(masked, cmap="YlOrRd", vmin=0, vmax=vmax, aspect="auto")
    plt.colorbar(im, ax=ax, shrink=0.8, label="Median |log₂FC|")
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, fill=True,
                                           color="#CCCCCC", zorder=2))
                ax.text(j, i, "—", ha="center", va="center",
                        fontsize=11, color="#666666", zorder=3)
            else:
                val = data_mat[i, j]
                tc  = "white" if val / vmax > 0.6 else "black"
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=8, color=tc, zorder=3)
    ax.set_title(f"{title}\n(n={n_genes:,} genes, min CPM ≥ {MIN_EXPR})", fontsize=10)
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=8)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(
    "Pairwise Median |log₂FC| — BAM CPM normalization\n"
    "(per million reads within BED regions)",
    fontsize=12, fontweight="bold", y=1.01
)

datasets = [
    ("Promoters (TSS ±500 bp)",   "promoter_bam_summary.npz"),
    ("Gene Body + 1 kb upstream", "genebody_bam_summary.npz"),
]

for ax, (name, npz) in zip(axes, datasets):
    cpm     = load_cpm(npz)
    lfc_mat = build_median_abslfc(cpm)
    n_genes = cpm.shape[0]
    vmax    = np.nanmax(np.where(np.eye(n, dtype=bool), np.nan, lfc_mat))
    draw(ax, lfc_mat, vmax, name, n_genes)

    print(f"\n{name} — Median |log₂FC|:")
    print(f"  {'':20}", "  ".join(f"{l:>13}" for l in labels))
    for i, li in enumerate(labels):
        row = "  ".join(f"{lfc_mat[i,j]:>13.3f}" if i != j else f"{'—':>13}" for j in range(n))
        print(f"  {li:<20}  {row}")

plt.tight_layout()
plt.savefig("bam_median_abslfc.svg", bbox_inches="tight")
print("\nSaved bam_median_abslfc.svg")
