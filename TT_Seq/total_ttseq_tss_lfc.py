"""
Log2 fold-change of total (unstranded) TT-seq signal around TSS.
Two panels: log2(+aux/CTL) and log2(+aux/−aux), same y-axis.
− strand gene profiles are flipped so all genes point right.
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
    "CTL":  "data/ttseq_bigwigs/INTS1_minus_unstranded.bw",
    "-aux": "data/ttseq_bigwigs/INTS11_minus_unstranded.bw",
    "+aux": "data/ttseq_bigwigs/INTS11_plus_unstranded.bw",
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

def extract_profiles(bw_path):
    bw = pyBigWig.open(bw_path)
    bw_chroms = set(bw.chroms().keys())
    mats = []
    for g in genes:
        chrom = g['chrom']
        c = chrom if chrom in bw_chroms else ("chr"+chrom if "chr"+chrom in bw_chroms else None)
        if c is None: continue
        clen = bw.chroms(c)
        w_lo = g['tss'] - HALF; w_hi = g['tss'] + HALF
        if w_lo < 0 or w_hi > clen: continue
        vals = bw.stats(c, w_lo, w_hi, type="mean", nBins=N_BINS)
        if vals is None or len(vals) != N_BINS: continue
        arr = np.array([v if v is not None else np.nan for v in vals], dtype=float)
        if g['strand'] == '-': arr = arr[::-1]
        mats.append(arr)
    bw.close()
    return np.array(mats)

means = {}
for cond, path in BW.items():
    print(f"Extracting {cond}...")
    mat = extract_profiles(path)
    means[cond] = np.nanmean(mat, axis=0)
    print(f"  n={mat.shape[0]}")

x = np.linspace(-HALF/1000, HALF/1000, N_BINS)

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
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("Distance from TSS (kb)")
    ax.set_ylabel(f"Total TT-seq {ylabel}")
    ax.set_xlim(-HALF/1000, HALF/1000)
    ax.set_ylim(*ylim)
    ax.set_title(ylabel, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    f"Total TT-seq log2 fold-change around TSS\n"
    f"genes ≥{MIN_LEN//1000} kb (n={len(genes)}), ±{HALF//1000} kb, {BIN_SIZE} bp bins",
    fontweight="bold")
fig.tight_layout()
fig.savefig("total_ttseq_tss_lfc.png", dpi=150, bbox_inches="tight")
fig.savefig("total_ttseq_tss_lfc.svg", bbox_inches="tight")
print("Saved total_ttseq_tss_lfc.{png,svg}")
