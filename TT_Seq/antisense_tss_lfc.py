"""
Log2 fold-change of antisense TT-seq signal around TSS.
Two panels: log2(+aux/CTL) and log2(+aux/−aux), same y-axis.
"""
import numpy as np, pyBigWig
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE_PATH = "genes_01.bed"
HALF      = 20_000
BIN_SIZE  = 100
N_BINS    = 2 * HALF // BIN_SIZE
MIN_LEN   = 40_000

BW = {
    "CTL":  {"fwd": "data/ttseq_bigwigs/INTS1_minus_fwd.bw",
              "rev": "data/ttseq_bigwigs/INTS1_minus_rev.bw"},
    "-aux": {"fwd": "data/ttseq_bigwigs/INTS11_minus_fwd.bw",
              "rev": "data/ttseq_bigwigs/INTS11_minus_rev.bw"},
    "+aux": {"fwd": "data/ttseq_bigwigs/INTS11_plus_fwd.bw",
              "rev": "data/ttseq_bigwigs/INTS11_plus_rev.bw"},
}

genes = []
for line in open(GENE_PATH):
    parts = line.strip().split('\t')
    if len(parts) < 6: continue
    chrom = parts[0]; start, end = int(parts[1]), int(parts[2])
    strand = parts[5]
    if end - start < MIN_LEN: continue
    tss = start if strand == '+' else end
    genes.append(dict(chrom=chrom, tss=tss, strand=strand))

print(f"Genes ≥{MIN_LEN//1000} kb: {len(genes)}")

def extract_profiles(bw_fwd, bw_rev):
    mats = []
    bw_chroms_fwd = set(bw_fwd.chroms().keys())
    bw_chroms_rev = set(bw_rev.chroms().keys())
    for g in genes:
        chrom = g['chrom']; strand = g['strand']; tss = g['tss']
        if strand == '+':
            bw = bw_rev; flip = False
            c  = chrom if chrom in bw_chroms_rev else ("chr"+chrom if "chr"+chrom in bw_chroms_rev else None)
        else:
            bw = bw_fwd; flip = True
            c  = chrom if chrom in bw_chroms_fwd else ("chr"+chrom if "chr"+chrom in bw_chroms_fwd else None)
        if c is None: continue
        clen = bw.chroms(c)
        w_lo = tss - HALF; w_hi = tss + HALF
        if w_lo < 0 or w_hi > clen: continue
        vals = bw.stats(c, w_lo, w_hi, type="mean", nBins=N_BINS)
        if vals is None or len(vals) != N_BINS: continue
        arr = np.array([v if v is not None else np.nan for v in vals], dtype=float)
        if flip: arr = arr[::-1]
        mats.append(arr)
    return np.array(mats)

means = {}
for cond, paths in BW.items():
    print(f"Extracting {cond}...")
    bwf = pyBigWig.open(paths["fwd"])
    bwr = pyBigWig.open(paths["rev"])
    mat = extract_profiles(bwf, bwr)
    bwf.close(); bwr.close()
    means[cond] = np.nanmean(mat, axis=0)
    print(f"  n={mat.shape[0]}")

x = np.linspace(-HALF/1000, HALF/1000, N_BINS)

# log2 FC — add small pseudocount relative to min nonzero signal
pseudo = min(np.nanmin(means[c][means[c] > 0]) for c in means) * 0.01

lfc_ctl  = np.log2((means["+aux"] + pseudo) / (means["CTL"]  + pseudo))
lfc_naux = np.log2((means["+aux"] + pseudo) / (means["-aux"] + pseudo))

# common y-axis
ymax = max(np.nanmax(np.abs(lfc_ctl)), np.nanmax(np.abs(lfc_naux))) * 1.05
ylim = (-ymax, ymax)

fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

panels = [
    (axes[0], lfc_ctl,  "#d73027", "log2(+aux / CTL)"),
    (axes[1], lfc_naux, "#d73027", "log2(+aux / −aux)"),
]
for ax, lfc, color, ylabel in panels:
    ax.plot(x, lfc, color=color, lw=1.8)
    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("Distance from TSS (kb)")
    ax.set_ylabel(f"Antisense TT-seq {ylabel}")
    ax.set_xlim(-HALF/1000, HALF/1000)
    ax.set_ylim(*ylim)
    ax.set_title(ylabel, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    f"Antisense TT-seq log2 fold-change around TSS\n"
    f"genes ≥{MIN_LEN//1000} kb (n={len(genes)}), ±{HALF//1000} kb, {BIN_SIZE} bp bins",
    fontweight="bold")
fig.tight_layout()
fig.savefig("antisense_tss_lfc.png", dpi=150, bbox_inches="tight")
fig.savefig("antisense_tss_lfc.svg", bbox_inches="tight")
print("Saved antisense_tss_lfc.{png,svg}")
