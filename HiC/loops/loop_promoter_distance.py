"""
Compare distance from CTCF loop anchors to nearest expressed gene promoter
for down vs stable loops. Uses raw +aux/ctl 2-fold classification with str<=20 filter.
"""
import bisect
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

LOOP_PATH    = "data/loops/CTCFLoops.txt"
EXPR_PATH    = "pc_parental_ints1_expressed.bed"
BIN          = 25000
STR_MAX      = 20
DOWN_THRESH  = -1.0
UP_THRESH    =  1.0

# --- load expressed TSSes (plain chrom names) ---
tss = {}  # chrom -> sorted array of TSS positions
for line in open(EXPR_PATH):
    c = line.rstrip().split("\t")
    ch = c[0]
    strand = c[5] if len(c) > 5 else "+"
    pos = int(c[1]) if strand == "+" else int(c[2])
    tss.setdefault(ch, []).append(pos)
for ch in tss:
    tss[ch] = np.array(sorted(tss[ch]))

def nearest_tss_dist(chrom, anchor_start):
    """Distance from anchor center to nearest expressed TSS."""
    ch = chrom[3:] if chrom.startswith("chr") else chrom
    arr = tss.get(ch)
    if arr is None or len(arr) == 0:
        return np.nan
    center = anchor_start + BIN // 2
    idx = bisect.bisect_left(arr, center)
    candidates = []
    if idx < len(arr):   candidates.append(abs(arr[idx] - center))
    if idx > 0:          candidates.append(abs(arr[idx-1] - center))
    return min(candidates)

# --- load loops and classify ---
down_dists, stable_dists = [], []  # per-anchor distances

for line in open(LOOP_PATH):
    c = line.rstrip().split("\t")
    chrom = c[0]; a1, a2 = int(c[1]), int(c[2])
    ctl, noaux, aux = float(c[3]), float(c[4]), float(c[5])
    if ctl <= 0 or aux <= 0: continue
    if ctl > STR_MAX: continue  # strength filter
    lfc = np.log2(aux / ctl)
    if lfc > UP_THRESH: continue  # skip up loops

    d1 = nearest_tss_dist(chrom, a1)
    d2 = nearest_tss_dist(chrom, a2)

    if lfc < DOWN_THRESH:
        if np.isfinite(d1): down_dists.append(d1)
        if np.isfinite(d2): down_dists.append(d2)
    else:
        if np.isfinite(d1): stable_dists.append(d1)
        if np.isfinite(d2): stable_dists.append(d2)

down_dists   = np.array(down_dists)
stable_dists = np.array(stable_dists)

# --- stats ---
u_stat, p_val = stats.mannwhitneyu(down_dists, stable_dists, alternative="less")
ks_stat, ks_p = stats.ks_2samp(down_dists, stable_dists)

print(f"Down anchors   n={len(down_dists):,}  median={np.median(down_dists)/1000:.1f} kb")
print(f"Stable anchors n={len(stable_dists):,}  median={np.median(stable_dists)/1000:.1f} kb")
print(f"Mann-Whitney (down < stable): U={u_stat:.0f}, p={p_val:.2e}")
print(f"KS test: D={ks_stat:.4f}, p={ks_p:.2e}")

for cutoff_kb in [25, 50, 100, 200]:
    cutoff = cutoff_kb * 1000
    frac_d = (down_dists   <= cutoff).mean() * 100
    frac_s = (stable_dists <= cutoff).mean() * 100
    print(f"  Within {cutoff_kb:3d} kb: down {frac_d:.1f}%  stable {frac_s:.1f}%")

# --- plot CDF ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: CDF up to 500 kb
ax = axes[0]
for arr, label, color in [
    (down_dists,   f"Down ({len(down_dists)//2} loops)", "#d73027"),
    (stable_dists, f"Stable ({len(stable_dists)//2} loops)", "#4575b4"),
]:
    xs = np.sort(arr)
    ys = np.arange(1, len(xs)+1) / len(xs)
    ax.plot(xs / 1000, ys, color=color, lw=2, label=label)
ax.set_xlabel("Distance to nearest expressed TSS (kb)")
ax.set_ylabel("Cumulative fraction of anchors")
ax.set_xlim(0, 500)
ax.set_ylim(0, 1)
ax.legend(frameon=False)
ax.set_title(f"p = {p_val:.1e} (Mann-Whitney)")
ax.spines[['top','right']].set_visible(False)

# Panel B: boxplot / violin of log10 distances
ax = axes[1]
data_kb = [down_dists/1000, stable_dists/1000]
parts = ax.violinplot(data_kb, positions=[0, 1], showmedians=True, showextrema=False)
colors = ["#d73027", "#4575b4"]
for pc, color in zip(parts['bodies'], colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
parts['cmedians'].set_color('black')
parts['cmedians'].set_linewidth(2)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"Down\n(n={len(down_dists):,})", f"Stable\n(n={len(stable_dists):,})"])
ax.set_ylabel("Distance to nearest expressed TSS (kb)")
ax.set_ylim(0, 1000)
ax.spines[['top','right']].set_visible(False)

fig.suptitle("CTCF loop anchors: distance to nearest expressed promoter")
fig.tight_layout()
fig.savefig("loop_promoter_distance.png", dpi=150, bbox_inches="tight")
fig.savefig("loop_promoter_distance.svg", bbox_inches="tight")
print("\nSaved loop_promoter_distance.{png,svg}")
