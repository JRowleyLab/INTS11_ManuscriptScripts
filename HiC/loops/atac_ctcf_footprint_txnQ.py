"""
ATAC-seq CTCF footprint at loop anchor motif sites, split into TT-seq
transcription-change quintiles (+aux vs ctl).

For each CTCF motif (MA0139.1) that falls inside a loop anchor bin (str<=20),
measure local TT-seq log2FC (+aux/ctl) at the motif center (±2.5 kb).
Split motifs into 5 quintiles, aggregate Tn5 cut sites ±300 bp per quintile,
and plot CTL vs +aux traces for each.
"""
import subprocess, os, tempfile, bisect, pysam
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyBigWig

LOOP_PATH  = "data/loops/CTCFLoops.txt"
PEAK_PATH  = "data/loops/GSM3715311_CTCF_CT_peaks.bed"
GENOME     = "data/annotation/hg38.fa"
MOTIF_DB   = "data/annotation/JASPAR_CORE_2014_vertebrates.meme"
BAM_CTL    = "data/atac/CTL.merged.nodup.sorted.bam"
BAM_AUX    = "data/atac/INTS11_AID.merged.nodup.sorted.bam"
BW_TXN     = "data/ttseq_bigwigs/log2_int11pvsint1m.bw"
FIMO_BIN   = "fimo"

BIN      = 25000
STR_MAX  = 20
HALF     = 1000
TXN_WIN  = 2500   # ±bp for TT-seq measurement at motif center
N_Q      = 5
FIMO_Q   = 0.05

# ── 1. collect all loop anchor bins (str<=20) ────────────────────────────────
anchor_bins = set()
for line in open(LOOP_PATH):
    c = line.rstrip().split("\t")
    if float(c[3]) > STR_MAX or float(c[3]) <= 0: continue
    chrom = c[0]
    anchor_bins.add((chrom, int(c[1])))
    anchor_bins.add((chrom, int(c[2])))
print(f"Loop anchor bins (str<=20): {len(anchor_bins)}")

def anchor_key(chrom, pos):
    return (chrom, (pos // BIN) * BIN)

# ── 2. run bedtools getfasta + FIMO (or reuse saved output) ──────────────────
fimo_tsv = "ctcf_fimo_output/fimo.tsv"
if not os.path.exists(fimo_tsv):
    os.makedirs("ctcf_fimo_output", exist_ok=True)
    fa_path = "ctcf_fimo_output/peaks.fa"
    print("Running bedtools getfasta...")
    subprocess.run(
        ["bedtools", "getfasta", "-fi", GENOME, "-bed", PEAK_PATH, "-fo", fa_path],
        check=True, capture_output=True
    )
    print("Running FIMO...")
    subprocess.run(
        [FIMO_BIN, "--motif", "MA0139.1",
         "--thresh", str(FIMO_Q), "--qv-thresh",
         "--oc", "ctcf_fimo_output", MOTIF_DB, fa_path],
        check=True, capture_output=True
    )
    print("FIMO done.")
else:
    print("Reusing existing FIMO output.")

# ── 3. parse FIMO → (chrom, center, strand) for motifs on loop anchors ───────
anchor_motifs = []   # (chrom, center, strand)
with open(fimo_tsv) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("motif_id") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 6: continue
        seq_name    = cols[2]
        motif_start = int(cols[3])   # 1-based within peak
        motif_stop  = int(cols[4])
        strand      = cols[5]
        chrom, coords = seq_name.rsplit(":", 1)
        peak_start  = int(coords.split("-")[0])   # 0-based
        center      = peak_start + (motif_start - 1 + motif_stop) // 2
        if anchor_key(chrom, center) in anchor_bins:
            anchor_motifs.append((chrom, center, strand))

print(f"CTCF motifs on loop anchors: {len(anchor_motifs)}")

# ── 4. assign TT-seq log2FC to each motif ────────────────────────────────────
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

txn_scores = np.array([txn_at(ch, cen) for ch, cen, _ in anchor_motifs])
bw.close()

valid = np.isfinite(txn_scores)
motifs_v = [m for m, ok in zip(anchor_motifs, valid) if ok]
scores_v = txn_scores[valid]
print(f"Motifs with valid TT-seq score: {len(motifs_v)}")

# ── 5. split into quintiles ───────────────────────────────────────────────────
order  = np.argsort(scores_v)
groups = np.array_split(order, N_Q)
quintile_motifs = [[motifs_v[i] for i in g] for g in groups]
for qi, (g, mots) in enumerate(zip(groups, quintile_motifs)):
    lo, hi = scores_v[g[0]], scores_v[g[-1]]
    print(f"  Q{qi+1}: n={len(mots)}  TT-seq log2FC [{lo:.3f}, {hi:.3f}]")

# ── 6. extract Tn5 cut sites ─────────────────────────────────────────────────
def extract_cuts(bam_path, motif_list, half=HALF):
    bam   = pysam.AlignmentFile(bam_path, "rb")
    total = bam.mapped
    width = 2 * half
    mat   = np.zeros((len(motif_list), width), dtype=np.float32)
    for i, (chrom, center, strand) in enumerate(motif_list):
        for read in bam.fetch(chrom, max(0, center - half - 5), center + half + 5):
            if read.is_unmapped or read.is_duplicate: continue
            cut = (read.reference_start + 4) if not read.is_reverse else (read.reference_end - 5 - 1)
            rel = cut - (center - half)
            if 0 <= rel < width:
                mat[i, rel] += 1
    bam.close()
    mat *= (1e6 / total)
    for i, (_, _, strand) in enumerate(motif_list):
        if strand == "-":
            mat[i] = mat[i, ::-1]
    return mat

def smooth(arr, w=10):
    return np.convolve(arr, np.ones(w) / w, mode="same")

print("\nExtracting cut sites per quintile...")
profiles_ctl, profiles_aux = [], []
for qi, mots in enumerate(quintile_motifs):
    print(f"  Q{qi+1} CTL...")
    profiles_ctl.append(extract_cuts(BAM_CTL, mots).mean(axis=0))
    print(f"  Q{qi+1} +aux...")
    profiles_aux.append(extract_cuts(BAM_AUX, mots).mean(axis=0))

# ── 7. plot ───────────────────────────────────────────────────────────────────
x      = np.arange(-HALF, HALF)
cmap   = plt.cm.RdYlBu_r
qcols  = [cmap(i / (N_Q - 1)) for i in range(N_Q)]
colors = {"CTL": "#4575b4", "+aux": "#d73027"}

fig, axes = plt.subplots(1, N_Q, figsize=(16, 4), sharey=True)
for qi, ax in enumerate(axes):
    n  = len(quintile_motifs[qi])
    lo = scores_v[groups[qi][0]]; hi = scores_v[groups[qi][-1]]
    for cond, profiles, col in [("CTL", profiles_ctl, colors["CTL"]),
                                 ("+aux", profiles_aux, colors["+aux"])]:
        ax.plot(x, smooth(profiles[qi]), color=col, lw=1.5, label=cond)
    ax.axvline(0, color="grey", lw=0.7, ls="--")
    ax.set_xlabel("Distance from CTCF motif center (bp)")
    if qi == 0:
        ax.set_ylabel("Mean Tn5 cut-site density (RPM)")
    ax.set_title(f"Q{qi+1}\nTT-seq [{lo:.2f}, {hi:.2f}]\nn={n}", fontsize=8)
    ax.legend(fontsize=7, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("ATAC-seq CTCF footprint by TT-seq quintile (log2 +aux/ctl)\nat loop anchors (str≤20)", fontweight="bold")
fig.tight_layout()
fig.savefig("atac_ctcf_footprint_txnQ.png", dpi=150, bbox_inches="tight")
fig.savefig("atac_ctcf_footprint_txnQ.svg", bbox_inches="tight")
print("\nSaved atac_ctcf_footprint_txnQ.{png,svg}")
