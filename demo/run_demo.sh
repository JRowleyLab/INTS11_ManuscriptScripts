#!/bin/bash
# run_demo.sh
#
# Runs both demos on the small chr19 dataset in demo/data/ and writes results
# to demo/output/. Run from the repository root:
#     bash demo/run_demo.sh
#
# Demo 1 (TSS heatmap) has two steps:
#   1a. stranded_computematrix.py : bigWig + gene BED -> per-gene TSS matrix (.npz)
#   1b. plot_ttseq_delta_heatmap.py : .npz -> heatmap + mean profile (.png/.svg)
# Demo 2 (scaled metaplot) is one step:
#   2.  Scaled1DMetaplotter_InfiniteAnchors.py : bigWig + anchor file -> scaled
#       average/median profiles, per-gene matrix and figures
#
# Each step is wrapped in `time` so the run time is printed.

set -e   # stop at the first command that fails

DATA=demo/data
OUT=demo/output
mkdir -p "$OUT"

echo "=== Demo 1a: build TSS-anchored matrix (20 kb up / 4 kb down, 50 bp bins) ==="
# The log2FC track is unstranded, so the same bigWig is given for both strands.
# --missing nan leaves bins without data out of the mean instead of counting them as 0.
time python HiC/preprocessing/ttseq/stranded_computematrix.py \
    --bed "$DATA/demo_genes.bed" \
    --plus-bw "$DATA/demo_log2FC_chr19.bw" \
    --minus-bw "$DATA/demo_log2FC_chr19.bw" \
    --upstream 20000 --downstream 4000 --bin-size 50 \
    --reference-point TSS --missing nan \
    --out "$OUT/demo_stranded.npz"

echo
echo "=== Demo 1b: plot heatmap + mean profile ==="
time python TT_Seq/plot_ttseq_delta_heatmap.py \
    --npz "$OUT/demo_stranded.npz" \
    --out-prefix "$OUT/demo_heatmap"

echo
echo "=== Demo 2: scaled metaplot (5 kb flank | gene body | 5 kb flank = 10,40,10 bins) ==="
# The output prefix ends in '_' because the script appends file suffixes directly.
time python Scaled1DMetaplotter_InfiniteAnchors.py \
    -i "$DATA/demo_log2FC_chr19.bw" \
    -a "$DATA/demo_anchors.txt" \
    -b 10,40,10 \
    -o "$OUT/demo_meta_" \
    -f -c darkred

echo
echo "Done. Outputs in $OUT/:"
ls -1 "$OUT"
