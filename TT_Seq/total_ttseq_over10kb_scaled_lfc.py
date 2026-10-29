"""
Log2 fold-change of total TT-seq around divergent gene pairs with >10 kb TSS gap.

Layout per gene:
  [left 10 kb flank | scaled inter-TSS region | right 10 kb flank]

The inter-TSS region is resampled to N_BODY=N_FLANK bins (equal visual width).
"""
import numpy as np, pyBigWig
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE_PATH = "over10kbgenes_01.bed"
FLANK     = 10_000
BIN_SIZE  = 100
N_FLANK   = FLANK // BIN_SIZE   # 100 bins per 10 kb

BW = {
    "CTL":  "data/ttseq_bigwigs/INTS1_minus_unstranded.bw",
    "-aux": "data/ttseq_bigwigs/INTS11_minus_unstranded.bw",
    "+aux": "data/ttseq_bigwigs/INTS11_plus_unstranded.bw",
}

regions = []
for line in open(GENE_PATH):
    parts = line.strip().split('\t')
    if len(parts) < 3: continue
    chrom = parts[0]; start, end = int(parts[1]), int(parts[2])
    gap = end - start
    if gap < 1: continue
    regions.append(dict(chrom=chrom, start=start, end=end, gap=gap))

mean_gap = np.mean([r['gap'] for r in regions])
N_BODY   = N_FLANK
N_TOTAL  = N_FLANK + N_BODY + N_FLANK

print(f"Regions: {len(regions)}")
print(f"Mean gap: {mean_gap/1000:.1f} kb → N_BODY={N_BODY} bins (equal to flanks, scaled)"
      f"  total: {N_TOTAL} bins")

def extract_profiles(bw_path):
    bw = pyBigWig.open(bw_path)
    bw_chroms = set(bw.chroms().keys())
    mats = []
    for r in regions:
        chrom = r['chrom']
        c = chrom if chrom in bw_chroms else ("chr"+chrom if "chr"+chrom in bw_chroms else None)
        if c is None: continue
        clen = bw.chroms(c)
        s, e = r['start'], r['end']
        if s - FLANK < 0 or e + FLANK > clen: continue

        def fetch(lo, hi, n):
            v = bw.stats(c, lo, hi, type="mean", nBins=n)
            if v is None or len(v) != n: return None
            return np.array([x if x is not None else np.nan for x in v], dtype=float)

        left  = fetch(s - FLANK, s, N_FLANK)
        body  = fetch(s, e,         N_BODY)
        right = fetch(e, e + FLANK, N_FLANK)
        if left is None or body is None or right is None: continue
        mats.append(np.concatenate([left, body, right]))
    bw.close()
    return np.array(mats)

means = {}
for cond, path in BW.items():
    print(f"Extracting {cond}...")
    mat = extract_profiles(path)
    means[cond] = np.nanmean(mat, axis=0)
    print(f"  n={mat.shape[0]}")

x       = np.arange(N_TOTAL, dtype=float)
body_kb = mean_gap / 1000

pseudo   = min(np.nanmin(means[c][means[c] > 0]) for c in means) * 0.01
lfc_ctl  = np.log2((means["+aux"] + pseudo) / (means["CTL"]  + pseudo))
lfc_naux = np.log2((means["+aux"] + pseudo) / (means["-aux"] + pseudo))

ymax = max(np.nanmax(np.abs(lfc_ctl)), np.nanmax(np.abs(lfc_naux))) * 1.05
ylim = (-ymax, ymax)

fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

for ax, lfc, ylabel in [
    (axes[0], lfc_ctl,  "log2(+aux / CTL)"),
    (axes[1], lfc_naux, "log2(+aux / −aux)"),
]:
    ax.plot(x, lfc, color="#d73027", lw=1.8)
    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.axvline(N_FLANK,          color="#555555", lw=1.0, ls=":")
    ax.axvline(N_FLANK + N_BODY, color="#555555", lw=1.0, ls=":")
    ax.axvspan(N_FLANK, N_FLANK + N_BODY, color="grey", alpha=0.10, zorder=0)
    ax.set_ylabel(f"Total TT-seq {ylabel}")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(*ylim)
    ax.set_title(ylabel, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ticks  = [0, N_FLANK, N_FLANK + N_BODY, N_TOTAL]
    labels = [f"−{FLANK//1000} kb", "TSS₁", f"TSS₂\n(scaled {body_kb:.0f} kb)", f"+{FLANK//1000} kb"]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)

fig.suptitle(
    f"Total TT-seq log2 fold-change — divergent pairs, inter-TSS gap >10 kb\n"
    f"(n={len(regions)} pairs, mean gap {body_kb:.0f} kb, "
    f"flanks ±{FLANK//1000} kb, {BIN_SIZE} bp bins)",
    fontweight="bold")
fig.tight_layout()
fig.savefig("total_ttseq_over10kb_scaled_lfc.png", dpi=150, bbox_inches="tight")
fig.savefig("total_ttseq_over10kb_scaled_lfc.svg", bbox_inches="tight")
print("Saved total_ttseq_over10kb_scaled_lfc.{png,svg}")
