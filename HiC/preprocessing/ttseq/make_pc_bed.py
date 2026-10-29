"""
Writes pc_ints1_expressed.bed  — protein-coding genes with INTS1_minus genebody CPM > 1.
Also writes upstream_distances.tsv with the gap (bp) from each gene's TSS
to the end of the nearest upstream gene on the same chromosome (any biotype).
"""
import numpy as np
import subprocess, re

GTF = 'data/annotation/gencode.v49.basic.annotation.gtf'
GENEBODY_NPZ = 'genebody_bam_summary.npz'
BED_GENEBODY = 'genebody_1kb_upstream.bed'

# ── 1. parse GTF: all genes (for upstream-distance calc) ─────────────────────
all_genes = {}   # gene_name -> (chrom, start0, end, strand)  [0-based half-open]
pc_genes  = set()

with open(GTF) as fh:
    for line in fh:
        if line.startswith('#') or '\tgene\t' not in line:
            continue
        m_type = re.search(r'gene_type "([^"]+)"', line)
        m_name = re.search(r'gene_name "([^"]+)"', line)
        if not m_name:
            continue
        gname = m_name.group(1)
        gtype = m_type.group(1) if m_type else ''
        cols  = line.split('\t')
        chrom = cols[0].lstrip('chr')
        start = int(cols[3]) - 1   # GTF is 1-based start → 0-based
        end   = int(cols[4])
        strand = cols[6]
        # keep last occurrence if duplicated (shouldn't matter much)
        all_genes[gname] = (chrom, start, end, strand)
        if gtype == 'protein_coding':
            pc_genes.add(gname)

print(f'All genes in GTF: {len(all_genes):,}')
print(f'Protein-coding: {len(pc_genes):,}')

# ── 2. load genebody BAM CPM, filter INTS1_minus > 1 ────────────────────────
INTS1_COL = 2   # PAR_minus=0, PAR_plus=1, INTS1_minus=2, INTS11_minus=3, INTS11_plus=4

d = np.load(GENEBODY_NPZ, allow_pickle=True)
mat = d['matrix']
valid_mask = ~np.isnan(mat).any(axis=1)
mat_v = mat[valid_mask]
totals = mat_v.sum(axis=0)
cpm    = mat_v / (totals / 1e6)

# BED gene names in same order as matrix rows (after NaN filter)
bed_names = []
with open(BED_GENEBODY) as f:
    for line in f:
        bed_names.append(line.split('\t')[3])
bed_names_valid = [g for g, v in zip(bed_names, valid_mask) if v]

ints1_cpm = cpm[:, INTS1_COL]
keep_mask = np.array([
    (g in pc_genes) and (ints1_cpm[i] > 1.0)
    for i, g in enumerate(bed_names_valid)
])
kept_names = [g for g, k in zip(bed_names_valid, keep_mask) if k]
print(f'Protein-coding + INTS1_minus > 1 CPM: {len(kept_names):,}')

# ── 3. write BED (original gene coords from GTF) ────────────────────────────
missing = 0
records = []
for gname in kept_names:
    if gname not in all_genes:
        missing += 1
        continue
    chrom, start, end, strand = all_genes[gname]
    records.append((chrom, start, end, gname, '.', strand))

with open('pc_ints1_expressed.bed', 'w') as out:
    for r in records:
        out.write('\t'.join(str(x) for x in r) + '\n')
print(f'Wrote pc_ints1_expressed.bed ({len(records):,} genes, {missing} missing from GTF)')

# ── 4. compute upstream gene distances ──────────────────────────────────────
# For each filtered gene, find the nearest upstream gene on same chrom.
# + strand gene TSS = start; upstream = any gene whose end <= start, same chrom
# - strand gene TSS = end;   upstream = any gene whose start >= end, same chrom
# Distance = gap between TSS and the upstream gene's closest boundary.

# Build per-chrom sorted structures: list of (end, start) for all genes
from collections import defaultdict
import bisect

chrom_plus_ends  = defaultdict(list)   # for searching upstream of + genes
chrom_minus_starts = defaultdict(list) # for searching upstream of - genes

for gname, (chrom, start, end, strand) in all_genes.items():
    chrom_plus_ends[chrom].append(end)         # end of any gene on this chrom
    chrom_minus_starts[chrom].append(start)    # start of any gene on this chrom

for chrom in chrom_plus_ends:
    chrom_plus_ends[chrom].sort()
for chrom in chrom_minus_starts:
    chrom_minus_starts[chrom].sort()

rows = []
for chrom, start, end, gname, dot, strand in records:
    if strand == '+':
        tss = start
        ends = chrom_plus_ends.get(chrom, [])
        # largest end strictly < tss
        idx = bisect.bisect_left(ends, tss) - 1
        if idx >= 0:
            dist = tss - ends[idx]
        else:
            dist = -1   # no upstream gene found
    else:
        tss = end
        starts = chrom_minus_starts.get(chrom, [])
        # smallest start strictly > tss
        idx = bisect.bisect_right(starts, tss)
        if idx < len(starts):
            dist = starts[idx] - tss
        else:
            dist = -1
    rows.append((gname, chrom, tss, strand, dist))

with open('upstream_distances.tsv', 'w') as out:
    out.write('gene\tchrom\ttss\tstrand\tupstream_dist_bp\n')
    for r in rows:
        out.write('\t'.join(str(x) for x in r) + '\n')
print('Wrote upstream_distances.tsv')
