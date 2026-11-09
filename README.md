# INTS11 manuscript analysis scripts

Analysis code for the INTS11 / Integrator manuscript (Nature Communications, *citation to be added*).

This README documents installation and a runnable demo for two tools:

| Tool | What it does | Manuscript figure |
|---|---|---|
| `TT_Seq/plot_ttseq_delta_heatmap.py` (with `HiC/preprocessing/ttseq/stranded_computematrix.py`) | Builds a strand-aware, TSS-anchored signal matrix from a bigWig and a gene list, then draws a per-gene heatmap with a mean ± SEM profile above it. | *Fig. X* |
| `Scaled1DMetaplotter_InfiniteAnchors.py` | Scaled 1D metaplot: splits each region into intervals between any number of anchors (e.g. upstream flank → TSS → TES → downstream flank), rescales each interval to a fixed number of bins, and averages across regions. | *Fig. X* |

The other scripts in this repository are the per-figure analysis code and are not covered by this demo.

---

## 1. System requirements

**Operating systems**
- Linux: tested on **Ubuntu 24.04.2 LTS** (x86-64).
- macOS: expected to work (all dependencies have macOS builds) but not tested.
- Windows: use **WSL2** (Ubuntu). `pyBigWig` does not build natively on Windows.

**Hardware**: a standard desktop or laptop. No GPU or other special hardware is needed. Peak memory is about 170 MB for the demo and well under 1 GB for the full manuscript gene set.

**Software dependencies**

| Package | Version tested |
|---|---|
| Python | 3.10.0 |
| numpy | 2.2.6 |
| matplotlib | 3.10.9 |
| pyBigWig | 0.3.29 |
| seaborn | 0.13.2 (pulls in pandas 2.3.3) |
| tqdm | 4.70.1 |

The exact versions are pinned in `requirements.txt`. Any Python ≥ 3.9 with recent versions of these packages should work.

## 2. Installation

```bash
git clone <repository URL>
cd INTS11_ManuscriptScripts

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If `pip` cannot install `pyBigWig` (it needs a C compiler plus the zlib and libcurl headers when no prebuilt wheel exists), use conda instead:

```bash
conda create -n ints11 -c conda-forge -c bioconda python=3.10 numpy matplotlib pybigwig seaborn tqdm
conda activate ints11
```

**Typical install time**: under 1 minute on a normal desktop with a broadband connection (15 s measured with pip and prebuilt wheels).
The conda route takes 2–5 minutes.

## 3. Demo

### Demo data (`demo/data/`, 4.8 MB)

The demo data is a **chromosome 19 subset of the real manuscript data**:

| File | Contents |
|---|---|
| `demo_log2FC_chr19.bw` | TT-seq log2 fold-change, INTS11 depletion vs INTS1 depletion, chr19 only. Base-level values are identical to the manuscript track. |
| `demo_genes.bed` | 623 expressed protein-coding genes on chr19 (BED6, strand in column 6), in the manuscript order: sorted by distance to the nearest upstream gene. |
| `demo_anchors.txt` | Anchor file for the metaplotter: `chrom, TSS−5kb, gene start, gene end, TES+5kb, strand`. 601 of the same genes (genes < 2 kb are excluded). |

These files were made from the full inputs with `demo/make_demo_subset.py`. That script is included to show where the files came from; you don't need to run it.

### Run the demo

From the repository root, with the environment activated:

```bash
bash demo/run_demo.sh
```

The script runs three commands, also listed here so you can run them by hand:

```bash
# 1a. Build the TSS-anchored matrix: 20 kb upstream, 4 kb downstream, 50 bp bins = 480 bins
python HiC/preprocessing/ttseq/stranded_computematrix.py \
    --bed demo/data/demo_genes.bed \
    --plus-bw demo/data/demo_log2FC_chr19.bw --minus-bw demo/data/demo_log2FC_chr19.bw \
    --upstream 20000 --downstream 4000 --bin-size 50 --reference-point TSS --missing nan \
    --out demo/output/demo_stranded.npz

# 1b. Heatmap + mean profile
python TT_Seq/plot_ttseq_delta_heatmap.py \
    --npz demo/output/demo_stranded.npz --out-prefix demo/output/demo_heatmap

# 2. Scaled metaplot: flank | gene body | flank = 10 | 40 | 10 bins
python Scaled1DMetaplotter_InfiniteAnchors.py \
    -i demo/data/demo_log2FC_chr19.bw -a demo/data/demo_anchors.txt \
    -b 10,40,10 -o demo/output/demo_meta_ -f -c darkred
```

### Expected output

The outputs are written to `demo/output/`. Reference copies are in `demo/expected_output/` for comparison.

| File | Description |
|---|---|
| `demo_stranded.npz` | Matrix of 623 genes × 480 bins (key `matrix`), plus gene names, chromosomes, strands and window settings. |
| `demo_heatmap.png` / `.svg` | Top: mean ± SEM log2FC profile. Bottom: one row per gene in input order. Vertical line at the TSS. Colour scale −0.4 to 0.4. |
| `demo_meta_ave_med.txt` | Two tab-separated rows, `averages` and `medians`, each with 60 values (one per scaled bin). |
| `demo_meta_mat.txt` | Per-gene scaled matrix of 601 genes × 60 bins (`numpy.loadtxt`-readable). Minus-strand genes are flipped to 5′→3′. |
| `demo_meta_ave.png`, `demo_meta_med.png` | Average and median scaled profile. |
| `demo_meta_heatmap.png` | Per-gene heatmap of the scaled matrix. |

Expected console output includes:

```
Done. Matrix shape: 623 genes × 480 bins
Matrix: 623 genes × 480 bins
Color limits: -0.4 to 0.4
Saved demo/output/demo_heatmap.svg / .png
```

In the expected heatmap profile, the mean log2FC rises over the ~5 kb upstream of the TSS and peaks just downstream of it (about +0.025).
The scaled metaplot average peaks near the TSS (bin 10).

To check numerically that your run matches the reference:

```bash
python -c "import numpy as np; \
a=np.load('demo/output/demo_stranded.npz')['matrix']; b=np.load('demo/expected_output/demo_stranded.npz')['matrix']; \
print('heatmap matrix matches:', np.allclose(a, b, equal_nan=True)); \
print('metaplot matrix matches:', np.allclose(np.loadtxt('demo/output/demo_meta_mat.txt'), np.loadtxt('demo/expected_output/demo_meta_mat.txt')))"
```

### Expected run time

These times were measured on Ubuntu 24.04. All steps are single-threaded, so a normal desktop should take about the same time.

| Step | Time |
|---|---|
| 1a. `stranded_computematrix.py` (623 genes) | ~40 s |
| 1b. `plot_ttseq_delta_heatmap.py` | ~2 s |
| 2. `Scaled1DMetaplotter_InfiniteAnchors.py` (601 genes) | ~7 s |
| **Total demo** | **~1 minute** (allow up to 2–3 min on an older laptop) |

The full manuscript gene set (8,874 genes) takes about 10 minutes for step 1a.

## 4. Running on your own data

### Input formats

- **bigWig**: any per-base signal (coverage, CPM, log2 ratio, etc.).
- **Gene BED** (for the heatmap): tab-separated BED6 (`chrom start end name score strand`) with the strand in column 6.
  The TSS is `start` for `+` genes and `end` for `−` genes. **Rows are plotted in file order**, so sort the BED the way you want the heatmap ordered.
  `stranded_computematrix.py` handles `chr`-prefix differences between the BED and the bigWig automatically.
- **Anchor file** (for the metaplotter): whitespace-separated `chrom anchor1 anchor2 … anchorN strand`.
  - Any N ≥ 2 is allowed, but every line must have the same number of anchors.
  - Coordinates must be **ascending on both strands**. The script flips `−` strand regions itself.
  - Chromosome names must match the bigWig **exactly** (`19` vs `chr19`).
  - Regions whose shortest interval is not longer than its bin count are skipped.

### Heatmap (two steps)

```bash
# Unstranded signal (e.g. a log2FC track): give the same bigWig twice
python HiC/preprocessing/ttseq/stranded_computematrix.py \
    --bed my_genes.bed --plus-bw my_signal.bw --minus-bw my_signal.bw \
    --upstream 20000 --downstream 4000 --bin-size 50 --reference-point TSS --missing nan \
    --out my_matrix.npz

# Stranded data: forward-strand bigWig for + genes, reverse for − genes (sense signal)
python HiC/preprocessing/ttseq/stranded_computematrix.py \
    --bed my_genes.bed --plus-bw my_fwd.bw --minus-bw my_rev.bw \
    --upstream 5000 --downstream 5000 --bin-size 25 --out my_matrix.npz

python TT_Seq/plot_ttseq_delta_heatmap.py --npz my_matrix.npz --out-prefix my_heatmap \
    --vmin -1 --vmax 1 --title-label "Treatment vs Control"
```

| `stranded_computematrix.py` option | Meaning (default) |
|---|---|
| `--bed`, `--plus-bw`, `--minus-bw`, `--out` | Required inputs and output `.npz`. |
| `--upstream`, `--downstream` | Window around the reference point in bp (2000 / 2000). |
| `--bin-size` | Bin size in bp (10). |
| `--reference-point` | `TSS`, `TES` or `center` (TSS). |
| `--missing` | Value for bins with no data; use `nan` to leave them out of the mean (0). |

| `plot_ttseq_delta_heatmap.py` option | Meaning (default) |
|---|---|
| `--npz` | Input matrix (`ttseq_int11pvsint1m_stranded.npz`). |
| `--out-prefix` | Writes `<prefix>.png` and `<prefix>.svg` (`ttseq_delta_INT11vsctl_heatmap`). |
| `--vmin`, `--vmax` | Colour limits (−0.4 / 0.4). |
| `--title-label` | Comparison name in the title and axis labels (`INTS11 vs INTS1`). |

The window and bin size are read from the `.npz`, and the x-axis ticks adjust to match.

### Scaled metaplot

```bash
python Scaled1DMetaplotter_InfiniteAnchors.py -i my_signal.bw -a my_anchors.txt \
    -b 10,40,10 -o my_meta_ -f -y 0,0.05 -c black -l 2
```

| Option | Meaning (default) |
|---|---|
| `-i` | Input bigWig (required). |
| `-a` | Anchor file (required). |
| `-b` | Bins per interval: one number for all intervals, or a comma-separated list with one value per interval (10). |
| `-o` | Output prefix. Suffixes are appended directly, so end it with `_` (`meta_out`). |
| `-f` | Also write PNG figures. Without it, only the text matrices are written. |
| `-y` | y-axis limits for the figures, `min,max` (auto). |
| `-c`, `-l` | Line colour and line width (black, 3). |

## 5. Reproducing the manuscript figures

The full-size inputs (`log2_int11pvsint1m.bw`, `genes_01_bydistance.bed`) are available from GEO (*accession to be added*).
`HiC/preprocessing/ttseq/run_stranded_matrices.sh` builds the matrices used in the manuscript.
`python TT_Seq/plot_ttseq_delta_heatmap.py` with no arguments then reproduces the published heatmap.

**Note on the demo subset**: the base-level values in `demo_log2FC_chr19.bw` are identical to the manuscript track.
bigWig files also contain precomputed summary levels, which pyBigWig uses for wide bins, and these are rebuilt when a subset is written.
As a result the demo heatmap matrix is bit-for-bit identical to running the full track on the same genes.
About 1% of metaplot bins differ by < 0.01 log2FC, and the average profile by < 1e-4.

## License and citation

*To be added.*
