#!/bin/bash
set -e

# Strand-aware TSS matrices (20 kb upstream / 4 kb downstream, 50 bp bins -> 480 bins, TSS at bin 400)
# sense: fwd bigwig on + genes, rev bigwig on - genes
# anti:  rev bigwig on + genes, fwd bigwig on - genes
# unstranded signals (Hi-C compartment delta, total TT-seq log2FC) pass the same bigwig twice

SCRIPT=preprocessing/ttseq/stranded_computematrix.py
BED=genes_01_bydistance.bed
ARGS="--bed $BED --upstream 20000 --downstream 4000 --bin-size 50 --reference-point TSS"
BW=data/ttseq_bigwigs

for cmp in int11pvsint11m int11pvsint1m int11pvsparm int11pvsparp
do
python3 $SCRIPT $ARGS --plus-bw $BW/log2_${cmp}_fwd.bw --minus-bw $BW/log2_${cmp}_rev.bw \
  --out log2_${cmp}_stranded.npz > log2_${cmp}_stranded.log &
python3 $SCRIPT $ARGS --plus-bw $BW/log2_${cmp}_rev.bw --minus-bw $BW/log2_${cmp}_fwd.bw \
  --out log2_${cmp}_anti_stranded.npz > log2_${cmp}_anti_stranded.log &
done

python3 $SCRIPT $ARGS --plus-bw data/compartments/delta_INT11vsctl.bw --minus-bw data/compartments/delta_INT11vsctl.bw \
  --out hic_delta_stranded.npz &
python3 $SCRIPT $ARGS --plus-bw data/compartments/delta_INT11vsnoaux.bw --minus-bw data/compartments/delta_INT11vsnoaux.bw \
  --out hic_delta_noaux_stranded.npz &

for cmp in int11pvsint1m int11pvsint11m
do
python3 $SCRIPT $ARGS --plus-bw $BW/log2_${cmp}.bw --minus-bw $BW/log2_${cmp}.bw --missing nan \
  --out ttseq_${cmp}_stranded.npz &
done

wait && echo "stranded matrices done"
