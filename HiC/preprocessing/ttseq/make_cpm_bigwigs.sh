#!/bin/bash
# Unstranded CPM bigwigs from MAPQ>=30, depth-subsampled TT-seq BAMs
# (parameters from the deeptools bamCoverage run logs)
for s in PAR_minus PAR_plus INTS1_minus INTS11_minus INTS11_plus
do
bamCoverage -b data/ttseq_bams/${s}_q30_sorted.bam -o data/ttseq_bigwigs/${s}_unstranded.bw \
  --normalizeUsing CPM --binSize 10 -p 8 &
done
wait
