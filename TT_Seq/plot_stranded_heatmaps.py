"""
Plot stranded TSS heatmaps + profile plots for all 4 log2FC comparisons.
Blue → white → orange, ±0.1, sorted by upstream gene distance.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

UPSTREAM   = 20000
DOWNSTREAM = 4000
BINSIZE    = 50
N_BINS     = (UPSTREAM + DOWNSTREAM) // BINSIZE   # 480
TSS_BIN    = UPSTREAM // BINSIZE                   # 400
VLIM       = 0.1

comparisons = [
    ("log2_int11pvsint11m_stranded.npz", "log2_int11pvsint11m_anti_stranded.npz", "INTS11+aux vs INTS11−aux"),
    ("log2_int11pvsint1m_stranded.npz",  "log2_int11pvsint1m_anti_stranded.npz",  "INTS11+aux vs INTS1−aux"),
    ("log2_int11pvsparm_stranded.npz",   "log2_int11pvsparm_anti_stranded.npz",   "INTS11+aux vs PAR−aux"),
    ("log2_int11pvsparp_stranded.npz",   "log2_int11pvsparp_anti_stranded.npz",   "INTS11+aux vs PAR+aux"),
]

cmap = mcolors.LinearSegmentedColormap.from_list(
    "blue_white_orange",
    ["#4393C3", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#E08214"],
    N=512,
)
SENSE_COLOR = "#E08214"
ANTI_COLOR  = "#4393C3"

# ── upstream distances ────────────────────────────────────────────────────────
dist_map = {}
with open("upstream_distances.tsv") as fh:
    next(fh)
    for line in fh:
        parts = line.strip().split("\t")
        dist_map[parts[0]] = int(parts[4])

# ── load matrices ─────────────────────────────────────────────────────────────
matrices = []
for sense_path, anti_path, label in comparisons:
    d      = np.load(sense_path, allow_pickle=True)
    mat_s  = d["matrix"]
    names  = d["names"]
    dists  = np.array([dist_map.get(n, -1) for n in names])
    d2     = np.load(anti_path, allow_pickle=True)
    mat_a  = d2["matrix"]
    matrices.append((mat_s, mat_a, dists, label))

n_genes = matrices[0][0].shape[0]

# ── pre-compute profiles and global y-limits ──────────────────────────────────
profiles = [(np.nanmean(ms, axis=0), np.nanmean(ma, axis=0))
            for ms, ma, _, _ in matrices]

all_sense = np.concatenate([ps for ps, _ in profiles])
all_anti  = np.concatenate([pa for _, pa in profiles])
sense_ylim = (0, 0.5)
anti_ylim  = (0, 0.5)

# ── x-axis tick positions (shared by profiles and heatmaps) ───────────────────
x_ticks  = [0, 100, 200, 300, TSS_BIN, TSS_BIN + 40, TSS_BIN + 80]
x_labels = ["-20kb", "-15kb", "-10kb", "-5kb", "TSS", "+2kb", "+4kb"]
bin_axis = np.arange(N_BINS)

# ── 50 % peak distance helper ─────────────────────────────────────────────────
def half_peak_distance(profile):
    """
    Find how far upstream of TSS the antisense signal extends before dropping
    to 50% of its upstream peak.

    The profile peaks at some upstream bin, then decays toward far-upstream
    baseline. We search leftward from the peak for the first bin where signal
    falls to <= 50% of peak, and report that bin's distance from the TSS.
    Returns (distance_bp, peak_bin, half_bin), or (None, None, None).
    """
    upstream = profile[:TSS_BIN]
    peak_bin = int(np.argmax(upstream))
    peak_val = upstream[peak_bin]
    half     = peak_val * 0.5

    # search leftward from peak for last bin (closest to peak) where signal <= half
    left_of_peak = upstream[:peak_bin]
    if left_of_peak.size == 0:
        return None, peak_bin, None
    candidates = np.where(left_of_peak <= half)[0]
    if len(candidates) == 0:
        return None, peak_bin, None

    half_bin    = int(candidates[-1])
    distance_bp = (TSS_BIN - half_bin) * BINSIZE
    return distance_bp, peak_bin, half_bin

# ── figure layout: 4 rows × 4 cols, profile rows short ───────────────────────
# rows: sense_profile, sense_heatmap, anti_profile, anti_heatmap
fig = plt.figure(figsize=(22, 22))
gs  = gridspec.GridSpec(
    4, 4,
    height_ratios=[1, 5, 1, 5],
    hspace=0.08, wspace=0.05,
)

print(f"\n{'Comparison':<30}  {'Anti 50%-peak distance':>25}")
print("-" * 57)

im_ref = None
for col, ((mat_s, mat_a, dists, label), (prof_s, prof_a)) in enumerate(zip(matrices, profiles)):

    # ── sense profile ─────────────────────────────────────────────────────────
    ax_sp = fig.add_subplot(gs[0, col])
    ax_sp.fill_between(bin_axis, prof_s, alpha=0.3, color=SENSE_COLOR)
    ax_sp.plot(bin_axis, prof_s, color=SENSE_COLOR, linewidth=0.8)
    ax_sp.axvline(TSS_BIN, color="black", linewidth=0.8, linestyle="--")
    ax_sp.axhline(0, color="grey", linewidth=0.4)
    ax_sp.set_xlim(0, N_BINS)
    ax_sp.set_ylim(*sense_ylim)
    ax_sp.set_xticks([])
    ax_sp.tick_params(labelsize=6)
    ax_sp.set_title(label, fontsize=8, fontweight="bold")
    if col == 0:
        ax_sp.set_ylabel("mean\nlog₂FC", fontsize=6)
        ax_sp.annotate("Sense", xy=(-0.22, 0.5), xycoords="axes fraction",
                        fontsize=9, fontweight="bold", rotation=90,
                        va="center", ha="center")
    else:
        ax_sp.set_yticks([])

    # ── sense heatmap ─────────────────────────────────────────────────────────
    ax_sh = fig.add_subplot(gs[1, col])
    im = ax_sh.imshow(mat_s, aspect="auto", cmap=cmap, vmin=-VLIM, vmax=VLIM,
                      interpolation="nearest", extent=[0, N_BINS, n_genes, 0])
    ax_sh.axvline(TSS_BIN, color="black", linewidth=0.8)
    xs = [TSS_BIN - d / BINSIZE for d in dists if 0 <= d <= UPSTREAM]
    ys = [i + 0.5 for i, d in enumerate(dists) if 0 <= d <= UPSTREAM]
    if xs:
        ax_sh.scatter(xs, ys, c="black", s=0.4, linewidths=0, zorder=5, rasterized=True)
    ax_sh.set_xticks(x_ticks)
    ax_sh.set_xticklabels(x_labels, fontsize=6, rotation=45, ha="right")
    ax_sh.set_xlim(0, N_BINS)
    ax_sh.set_xlabel("Distance from TSS", fontsize=6)
    ax_sh.tick_params(labelsize=6)
    if col == 0:
        ax_sh.set_ylabel("Genes", fontsize=7)
    else:
        ax_sh.set_yticks([])
    im_ref = im

    # ── antisense profile ─────────────────────────────────────────────────────
    ax_ap = fig.add_subplot(gs[2, col])
    ax_ap.fill_between(bin_axis, prof_a, alpha=0.3, color=ANTI_COLOR)
    ax_ap.plot(bin_axis, prof_a, color=ANTI_COLOR, linewidth=0.8)
    ax_ap.axvline(TSS_BIN, color="black", linewidth=0.8, linestyle="--")
    ax_ap.axhline(0, color="grey", linewidth=0.4)
    ax_ap.set_xlim(0, N_BINS)
    ax_ap.set_ylim(*anti_ylim)
    ax_ap.set_xticks([])
    ax_ap.tick_params(labelsize=6)
    if col == 0:
        ax_ap.set_ylabel("mean\nlog₂FC", fontsize=6)
        ax_ap.annotate("Antisense", xy=(-0.22, 0.5), xycoords="axes fraction",
                        fontsize=9, fontweight="bold", rotation=90,
                        va="center", ha="center")
    else:
        ax_ap.set_yticks([])

    # 50 % peak line on antisense profile
    d50, peak_bin, half_bin = half_peak_distance(prof_a)
    if d50 is not None:
        ax_ap.axvline(half_bin, color=ANTI_COLOR, linewidth=1.0,
                      linestyle=":", alpha=0.8)
        ax_ap.axvline(peak_bin, color=ANTI_COLOR, linewidth=0.8,
                      linestyle="-", alpha=0.5)
        ymax = ax_ap.get_ylim()[1]
        ax_ap.text(half_bin - 1, ymax * 0.9,
                   f"−{d50/1000:.0f}kb", fontsize=6, color=ANTI_COLOR,
                   ha="right", va="top")
        print(f"{label:<30}  {d50:>10,} bp  ({d50/1000:.1f} kb) upstream of TSS")
    else:
        print(f"{label:<30}  {'never reaches 50% of peak':>25}")

    # ── antisense heatmap ─────────────────────────────────────────────────────
    ax_ah = fig.add_subplot(gs[3, col])
    ax_ah.imshow(mat_a, aspect="auto", cmap=cmap, vmin=-VLIM, vmax=VLIM,
                 interpolation="nearest", extent=[0, N_BINS, n_genes, 0])
    ax_ah.axvline(TSS_BIN, color="black", linewidth=0.8)
    if xs:
        ax_ah.scatter(xs, ys, c="black", s=0.4, linewidths=0, zorder=5, rasterized=True)
    ax_ah.set_xticks(x_ticks)
    ax_ah.set_xticklabels(x_labels, fontsize=6, rotation=45, ha="right")
    ax_ah.set_xlim(0, N_BINS)
    ax_ah.set_xlabel("Distance from TSS", fontsize=6)
    ax_ah.tick_params(labelsize=6)
    if col == 0:
        ax_ah.set_ylabel("Genes", fontsize=7)
    else:
        ax_ah.set_yticks([])

fig.suptitle(
    "Stranded log₂FC around TSS  |  sorted by upstream gene distance",
    fontsize=11, fontweight="bold", y=1.005,
)

cbar = fig.colorbar(
    im_ref,
    ax=fig.axes,
    shrink=0.25, pad=0.01, fraction=0.01,
    label="log₂FC",
)
cbar.ax.tick_params(labelsize=7)

plt.savefig("stranded_heatmaps.svg", bbox_inches="tight")
plt.savefig("stranded_heatmaps.png", bbox_inches="tight", dpi=150)
print("\nSaved stranded_heatmaps.svg / .png")
