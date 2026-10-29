#!/usr/bin/env python3
"""
stranded_computematrix.py

Strand-aware alternative to deepTools computeMatrix reference-point.

For each region in the BED file:
  • + strand → coverage pulled from --plus-bw; array left = upstream, right = downstream
  • − strand → coverage pulled from --minus-bw; array is reversed so left = upstream,
                right = downstream (i.e., higher genomic coords come first)

The reference point (TSS, TES, or center) is also strand-aware.

Output: compressed .npz with keys
  matrix          float32 [n_genes × n_bins]
  names           str     [n_genes]
  chroms          str     [n_genes]
  strands         str     [n_genes]   '+' or '-'
  upstream        int32
  downstream      int32
  bin_size        int32
  reference_point str

Usage example
-------------
python stranded_computematrix.py \\
    --bed genes.bed \\
    --plus-bw sample_fwd.bw \\
    --minus-bw sample_rev.bw \\
    --upstream 5000 --downstream 2000 --bin-size 10 \\
    --reference-point TSS \\
    --out matrix.npz
"""

import argparse
import sys
import numpy as np
import pyBigWig


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bed',              required=True,
                   help='BED file; strand must be in column 6 (1-based)')
    p.add_argument('--plus-bw',          required=True,
                   help='BigWig for the forward / + strand')
    p.add_argument('--minus-bw',         required=True,
                   help='BigWig for the reverse / − strand')
    p.add_argument('--upstream',         type=int, default=2000,
                   help='bp upstream of reference point (default: 2000)')
    p.add_argument('--downstream',       type=int, default=2000,
                   help='bp downstream of reference point (default: 2000)')
    p.add_argument('--bin-size',         type=int, default=10,
                   help='Bin size in bp (default: 10)')
    p.add_argument('--reference-point',  choices=['TSS', 'TES', 'center'],
                   default='TSS',
                   help='Anchor point; strand-aware (default: TSS)')
    p.add_argument('--missing',          type=float, default=0.0,
                   help='Fill value for bins with no bigWig data (default: 0)')
    p.add_argument('--out',              required=True,
                   help='Output .npz file')
    return p.parse_args()


# ── BED parsing ────────────────────────────────────────────────────────────────

def parse_bed(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols   = line.rstrip('\n').split('\t')
            chrom  = cols[0]
            start  = int(cols[1])
            end    = int(cols[2])
            name   = cols[3] if len(cols) > 3 else f'{chrom}:{start}-{end}'
            strand = cols[5] if len(cols) > 5 else '+'
            if strand not in ('+', '-'):
                strand = '+'
            genes.append((chrom, start, end, name, strand))
    return genes


# ── Reference-point calculation ────────────────────────────────────────────────

def ref_coord(start, end, strand, ref_point):
    """
    Return the genomic anchor position (0-based).

    BED convention: start is 0-based inclusive, end is 0-based exclusive.
    For minus-strand genes the TSS is at the rightmost coordinate (= end in BED,
    matching deepTools behaviour).
    """
    if ref_point == 'TSS':
        return start if strand == '+' else end
    elif ref_point == 'TES':
        return end if strand == '+' else start
    else:  # center
        return (start + end) // 2


# ── BigWig fetch ───────────────────────────────────────────────────────────────

def _resolve_chrom(bw, chrom):
    """Return the chromosome name as it appears in the bigWig, or None."""
    if bw.chroms(chrom) is not None:
        return chrom
    alt = ("chr" + chrom) if not chrom.startswith("chr") else chrom[3:]
    return alt if bw.chroms(alt) is not None else None


def fetch_bins(bw, chrom, reg_start, reg_end, n_bins, missing_val):
    """
    Return a float32 array of length n_bins covering [reg_start, reg_end).

    Coordinates are clamped to [0, chrom_len); bins outside that range are
    filled with missing_val.  Bins where the bigWig has no data are also
    filled with missing_val.  Handles chr-prefix mismatches automatically.
    """
    chrom = _resolve_chrom(bw, chrom)
    out = np.full(n_bins, missing_val, dtype=np.float32)

    if chrom is None or reg_start >= reg_end:
        return out

    chrom_len = bw.chroms(chrom)

    span = reg_end - reg_start
    q0   = max(0, reg_start)
    q1   = min(chrom_len, reg_end)

    if q0 >= q1:
        return out

    # Map clamped genomic coords to bin indices (proportional)
    i0 = int(round((q0 - reg_start) / span * n_bins))
    i1 = int(round((q1 - reg_start) / span * n_bins))
    i0 = max(0, min(i0, n_bins))
    i1 = max(0, min(i1, n_bins))

    if i0 >= i1:
        return out

    n_query = i1 - i0
    vals = bw.stats(chrom, q0, q1, type='mean', nBins=n_query)
    if vals:
        for k, v in enumerate(vals):
            if v is not None:
                out[i0 + k] = float(v)

    return out


# ── Per-gene processing ────────────────────────────────────────────────────────

def process_gene(gene, plus_bw, minus_bw,
                 upstream, downstream, n_bins, missing_val, ref_point):
    chrom, start, end, name, strand = gene
    ref = ref_coord(start, end, strand, ref_point)

    if strand == '+':
        # Region: [ref − upstream,  ref + downstream)
        # Column 0 = upstream end, column -1 = downstream end — no flip needed.
        arr = fetch_bins(plus_bw, chrom,
                         ref - upstream, ref + downstream,
                         n_bins, missing_val)
    else:
        # For a − strand gene, "upstream" is at HIGHER genomic coords.
        # Region to fetch: [ref − downstream,  ref + upstream)
        # After fetching, arr[0] = lowest coord = downstream end.
        # Reverse so arr[0] = upstream end (5'→3' orientation).
        arr = fetch_bins(minus_bw, chrom,
                         ref - downstream, ref + upstream,
                         n_bins, missing_val)
        arr = arr[::-1].copy()

    return arr


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if (args.upstream + args.downstream) % args.bin_size != 0:
        print(f"Warning: (upstream + downstream) = {args.upstream + args.downstream} "
              f"is not evenly divisible by bin_size = {args.bin_size}. "
              f"Using floor division.", file=sys.stderr)

    n_bins = (args.upstream + args.downstream) // args.bin_size
    ref_bin = args.upstream // args.bin_size   # column index of the reference point

    print(f"Reading BED: {args.bed}")
    genes = parse_bed(args.bed)
    n_genes = len(genes)
    print(f"  {n_genes:,} regions")
    print(f"  reference-point : {args.reference_point}")
    print(f"  upstream        : {args.upstream:,} bp")
    print(f"  downstream      : {args.downstream:,} bp")
    print(f"  bin size        : {args.bin_size} bp")
    print(f"  total bins      : {n_bins}  (reference at bin {ref_bin})")

    plus_bw  = pyBigWig.open(args.plus_bw)
    minus_bw = pyBigWig.open(args.minus_bw)

    matrix  = np.full((n_genes, n_bins), args.missing, dtype=np.float32)
    names   = []
    chroms  = []
    strands = []

    print(f"Processing…", flush=True)
    dot_every = max(1, n_genes // 50)

    for i, gene in enumerate(genes):
        matrix[i] = process_gene(gene, plus_bw, minus_bw,
                                  args.upstream, args.downstream,
                                  n_bins, args.missing, args.reference_point)
        names.append(gene[3])
        chroms.append(gene[0])
        strands.append(gene[4])

        if (i + 1) % dot_every == 0 or (i + 1) == n_genes:
            pct = (i + 1) / n_genes * 100
            print(f"  {i+1:>{len(str(n_genes))}}/{n_genes}  ({pct:5.1f}%)", end='\r', flush=True)

    print()  # newline after progress

    plus_bw.close()
    minus_bw.close()

    print(f"Saving → {args.out}")
    np.savez_compressed(
        args.out,
        matrix          = matrix,
        names           = np.array(names,   dtype=object),
        chroms          = np.array(chroms,  dtype=object),
        strands         = np.array(strands, dtype=object),
        upstream        = np.int32(args.upstream),
        downstream      = np.int32(args.downstream),
        bin_size        = np.int32(args.bin_size),
        reference_point = np.str_(args.reference_point),
    )
    print(f"Done. Matrix shape: {matrix.shape[0]:,} genes × {matrix.shape[1]:,} bins")


if __name__ == '__main__':
    main()
