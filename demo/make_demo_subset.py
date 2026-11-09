#!/usr/bin/env python3
"""
make_demo_subset.py

One-off script that built the small demo dataset in demo/data/ from the full
manuscript inputs. It is kept in the repository so the provenance of the demo
files is transparent; you do NOT need to run it to use the demo.

It takes:
  - the genome-wide TT-seq log2 fold-change bigWig (INTS11 vs INTS1)
  - the ordered gene list used for the manuscript heatmaps (genes_01_bydistance.bed)
and keeps only one chromosome (chr19 by default), writing:
  - demo_log2FC_chr<CHROM>.bw : bigWig with just that chromosome
  - demo_genes.bed            : BED6 genes on that chromosome, in the original order
  - demo_anchors.txt          : 4-anchor file for Scaled1DMetaplotter_InfiniteAnchors.py
                                (TSS-5kb | gene body | TES+5kb), same genes

Note on bigWig zoom levels: the base-level values in the demo bigWig are
identical to the original. bigWig files also store pre-computed summaries
("zoom levels") that pyBigWig's stats() uses for wide bins, and these are
rebuilt when a new file is written. As a result the TSS heatmap matrix
(50 bp bins) is bit-for-bit identical to the full-file run, while about 1% of
metaplot bins differ by < 0.01 log2FC and the average profile by < 1e-4.

Usage (from the repository root):
    python demo/make_demo_subset.py \
        --bw TT_Seq/log2_int11pvsint1m.bw \
        --bed TT_Seq/genes_01_bydistance.bed \
        --chrom 19 --outdir demo/data
"""
import argparse
import os
import pyBigWig

# Flank added on each side of the gene body for the metaplot anchors (bp).
FLANK = 5000
# Genes shorter than this are left out of the anchor file so that every
# interval is much longer than the number of bins it is split into.
MIN_GENE_LEN = 2000


def parse_args():
    """Command-line options: input files, chromosome to keep, output folder."""
    p = argparse.ArgumentParser(description="Cut the manuscript inputs down to one chromosome for the demo.")
    p.add_argument("--bw", required=True, help="Genome-wide log2FC bigWig")
    p.add_argument("--bed", required=True, help="Ordered gene BED6 (genes_01_bydistance.bed)")
    p.add_argument("--chrom", default="19", help="Chromosome to keep, named as in the bigWig (default: 19)")
    p.add_argument("--outdir", default="demo/data", help="Output folder (default: demo/data)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    chrom = args.chrom

    # ── 1. bigWig: copy every interval on the chosen chromosome ────────────────
    src = pyBigWig.open(args.bw)
    chrom_len = src.chroms(chrom)
    intervals = src.intervals(chrom)   # tuple of (start, end, value)
    src.close()

    bw_out_path = os.path.join(args.outdir, f"demo_log2FC_chr{chrom}.bw")
    out = pyBigWig.open(bw_out_path, "w")
    # The header lists only the one chromosome we keep.
    out.addHeader([(chrom, chrom_len)])
    # addEntries takes parallel lists of chrom names, starts, ends and values.
    out.addEntries([chrom] * len(intervals),
                   [iv[0] for iv in intervals],
                   ends=[iv[1] for iv in intervals],
                   values=[float(iv[2]) for iv in intervals])
    out.close()

    # ── 2. BED: keep genes on this chromosome, preserving the original order ───
    # The original file is sorted by distance to the nearest upstream gene,
    # and the heatmap plots rows in file order, so we must not re-sort.
    genes = []
    with open(args.bed) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if cols[0] == chrom:
                genes.append(cols)

    bed_out_path = os.path.join(args.outdir, "demo_genes.bed")
    with open(bed_out_path, "w") as fh:
        for cols in genes:
            fh.write("\t".join(cols) + "\n")

    # ── 3. Anchor file for the scaled metaplotter ──────────────────────────────
    # Format: chrom  a1  a2  a3  a4  strand, coordinates ascending on BOTH strands
    # (the metaplotter flips minus-strand genes itself).
    # a1 = start - FLANK, a2 = start, a3 = end, a4 = end + FLANK.
    anchors_out_path = os.path.join(args.outdir, "demo_anchors.txt")
    n_anchors = 0
    with open(anchors_out_path, "w") as fh:
        for cols in genes:
            start, end, strand = int(cols[1]), int(cols[2]), cols[5]
            # skip short genes and genes whose flanks would run off the chromosome
            if end - start < MIN_GENE_LEN or start - FLANK < 0 or end + FLANK > chrom_len:
                continue
            fh.write(f"{chrom}\t{start - FLANK}\t{start}\t{end}\t{end + FLANK}\t{strand}\n")
            n_anchors += 1

    # ── 4. Report what was written (numbers go into the README) ────────────────
    print(f"chromosome {chrom}: {chrom_len:,} bp, {len(intervals):,} bigWig intervals")
    print(f"{bw_out_path}: {os.path.getsize(bw_out_path) / 1e6:.1f} MB")
    print(f"{bed_out_path}: {len(genes)} genes")
    print(f"{anchors_out_path}: {n_anchors} genes (>= {MIN_GENE_LEN} bp, flanks inside chromosome)")


if __name__ == "__main__":
    main()
