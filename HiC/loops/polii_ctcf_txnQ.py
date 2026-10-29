"""
PolII Ser2P and Ser5P ChIP-seq signal across CTCF motif sites at loop anchors,
split into TT-seq quintiles (log2 +aux/ctl). Mirrors atac_ctcf_footprint_txnQ.py.
Two rows (Ser2, Ser5), five columns (Q1-Q5), CTL vs +aux per panel.
Y-scale is uniform across all panels.
"""
import os, bisect
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyBigWig

LOOP_PATH = "data/loops/CTCFLoops.txt"
FIMO_TSV  = "ctcf_fimo_output/fimo.tsv"
BW_TXN    = "data/ttseq_bigwigs/log2_int11pvsint1m.bw"
S2P_DIR   = "data/cuttag"

MARKS = [
    ("Ser2P", f"{S2P_DIR}/S2P_INTS1-aux_merged.bw",  f"{S2P_DIR}/S2P_INTS11+aux_merged.bw"),
    ("Ser5P", f"{S2P_DIR}/S5P_INTS1-aux_merged.bw",  f"{S2P_DIR}/S5P_INTS11+aux_merged.bw"),
]

BIN     = 25000
STR_MAX = 20
HALF    = 1000
TXN_WIN = 2500
N_Q     = 5
BIN_BP  = 10   # resolution for bigwig extraction

# ── 1. anchor bins ────────────────────────────────────────────────────────────
anchor_bins = set()
for line in open(LOOP_PATH):
    c = line.rstrip().split("\t")
    if float(c[3]) > STR_MAX or float(c[3]) <= 0: continue
    chrom = c[0]
    anchor_bins.add((chrom, int(c[1])))
    anchor_bins.add((chrom, int(c[2])))

def anchor_key(chrom, pos):
    return (chrom, (pos // BIN) * BIN)

# ── 2. parse FIMO → motifs on anchors ────────────────────────────────────────
anchor_motifs = []
with open(FIMO_TSV) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("motif_id") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 6: continue
        seq_name    = cols[2]
        motif_start = int(cols[3])
        motif_stop  = int(cols[4])
        strand      = cols[5]
        chrom, coords = seq_name.rsplit(":", 1)
        peak_start  = int(coords.split("-")[0])
        center      = peak_start + (motif_start - 1 + motif_stop) // 2
        if anchor_key(chrom, center) in anchor_bins:
            anchor_motifs.append((chrom, center, strand))

print(f"CTCF motifs on loop anchors: {len(anchor_motifs)}")

# ── 3. TT-seq scores → quintiles ─────────────────────────────────────────────
bw = pyBigWig.open(BW_TXN)
bw_chroms = set(bw.chroms().keys())

def txn_at(chrom, center):
    for c in (chrom, chrom.replace("chr","") if chrom.startswith("chr") else "chr"+chrom):
        if c in bw_chroms:
            clen = bw.chroms(c)
            s, e = max(0, center - TXN_WIN), min(clen, center + TXN_WIN)
            v = bw.stats(c, s, e, type="mean")
            return float(v[0]) if v and v[0] is not None else np.nan
    return np.nan

scores = np.array([txn_at(ch, cen) for ch, cen, _ in anchor_motifs])
bw.close()

valid     = np.isfinite(scores)
motifs_v  = [m for m, ok in zip(anchor_motifs, valid) if ok]
scores_v  = scores[valid]
order     = np.argsort(scores_v)
groups    = np.array_split(order, N_Q)
quintile_motifs = [[motifs_v[i] for i in g] for g in groups]
for qi, (g, mots) in enumerate(zip(groups, quintile_motifs)):
    lo, hi = scores_v[g[0]], scores_v[g[-1]]
    print(f"  Q{qi+1}: n={len(mots)}  TT-seq [{lo:.3f}, {hi:.3f}]")

# ── 4. extract bigwig profile around motif centers ───────────────────────────
def bw_profile(bw_path, motif_list, half=HALF, bin_bp=BIN_BP):
    """
    Return (n_motifs, 2*half) array of mean signal per bp position.
    Extracts in bin_bp-wide bins then upsamples for consistency with footprint plots.
    Reverse-strand motifs are flipped.
    """
    bw   = pyBigWig.open(bw_path)
    bwcs = set(bw.chroms().keys())
    width = 2 * half
    n_bins = width // bin_bp
    mat  = np.zeros((len(motif_list), n_bins), dtype=np.float32)

    for i, (chrom, center, strand) in enumerate(motif_list):
        c = chrom if chrom in bwcs else (chrom.replace("chr","") if chrom.startswith("chr") else "chr"+chrom)
        if c not in bwcs: continue
        clen = bw.chroms(c)
        s    = max(0, center - half)
        e    = min(clen, center + half)
        if e <= s: continue
        vals = bw.stats(c, s, e, type="mean", nBins=n_bins)
        if vals is None: continue
        arr = np.array([v if v is not None else 0.0 for v in vals], dtype=np.float32)
        # pad if region was clipped at chrom boundary
        if len(arr) < n_bins:
            full = np.zeros(n_bins, dtype=np.float32)
            offset = (s - (center - half)) // bin_bp
            full[offset:offset+len(arr)] = arr
            arr = full
        mat[i] = arr
        if strand == "-":
            mat[i] = mat[i, ::-1]

    bw.close()
    # upsample bins → per-bp for consistent x-axis with footprint plots
    return np.repeat(mat, bin_bp, axis=1)

def smooth(arr, w=10):
    return np.convolve(arr, np.ones(w) / w, mode="same")

# compute all profiles
all_profiles = {}   # (mark_label, cond) -> list of per-quintile mean profiles
for mark_label, bw_ctl, bw_aux in MARKS:
    for cond, bw_path in [("CTL", bw_ctl), ("+aux", bw_aux)]:
        print(f"  {mark_label} {cond}...")
        profs = []
        for mots in quintile_motifs:
            mat = bw_profile(bw_path, mots)
            profs.append(mat.mean(axis=0))
        all_profiles[(mark_label, cond)] = profs

# ── 5. plot ───────────────────────────────────────────────────────────────────
x      = np.arange(-HALF, HALF)
colors = {"CTL": "#4575b4", "+aux": "#d73027"}

ymin, ymax = 0, 2

fig, axes = plt.subplots(2, N_Q, figsize=(16, 7), sharey=True)

for row, (mark_label, _, _) in enumerate(MARKS):
    for qi in range(N_Q):
        ax  = axes[row, qi]
        n   = len(quintile_motifs[qi])
        lo  = scores_v[groups[qi][0]]; hi = scores_v[groups[qi][-1]]
        for cond in ("CTL", "+aux"):
            prof = all_profiles[(mark_label, cond)][qi]
            ax.plot(x, smooth(prof), color=colors[cond], lw=1.5, label=cond)
        ax.axvline(0, color="grey", lw=0.7, ls="--")
        ax.set_ylim(ymin, ymax)
        if qi == 0:
            ax.set_ylabel(f"{mark_label}\nMean ChIP signal")
        if row == 0:
            ax.set_title(f"Q{qi+1}\n[{lo:.2f}, {hi:.2f}]\nn={n}", fontsize=8)
        if row == len(MARKS) - 1:
            ax.set_xlabel("Distance from CTCF motif (bp)", fontsize=8)
        if row == 0 and qi == N_Q - 1:
            ax.legend(fontsize=7, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("PolII Ser2P and Ser5P at CTCF motifs on loop anchors\nby TT-seq quintile (log2 +aux/ctl)",
             fontweight="bold")
fig.tight_layout()
fig.savefig("polii_ctcf_txnQ.png", dpi=150, bbox_inches="tight")
fig.savefig("polii_ctcf_txnQ.svg", bbox_inches="tight")
print("\nSaved polii_ctcf_txnQ.{png,svg}")
