"""
Classify CTCF loops as down / stable / up from raw SIP loop strengths.

Loop change = log2(+aux / ctl) from CTCFLoops.txt (chr, a1, a2, ctl, -aux, +aux; 25 kb bins).
  down   : log2FC < -1
  up     : log2FC > +1
  stable : otherwise
Loops with ctl APA strength > 20 (ctl_totloops_strength.tab) are removed first
(pericentromeric / segmental-duplication artifacts).

No MA normalization: the -aux INTS11-AID condition is leaky, so the global loop
decrease in +aux vs ctl is treated as real.

Writes ctcf_loops_{down,stable,up}_raw2fold.bedpe
"""
import numpy as np

LOOPS    = "data/loops/CTCFLoops.txt"
STRENGTH = "data/loops/ctl_totloops_strength.tab"
RES      = 25000
MAX_STR  = 20

strength = {}
for line in open(STRENGTH):
    p = line.split()
    try:
        strength[(p[1], int(p[2]), int(p[5]))] = float(p[7])
    except (IndexError, ValueError):
        continue   # header

out = {"down": [], "stable": [], "up": []}
for line in open(LOOPS):
    p = line.split()
    chrom, a1, a2 = p[0], int(p[1]), int(p[2])
    ctl, noaux, aux = map(float, p[3:6])
    if strength.get((chrom, a1, a2), 0) > MAX_STR:
        continue
    lfc = np.log2(aux / ctl)
    cls = "down" if lfc < -1 else "up" if lfc > 1 else "stable"
    out[cls].append(f"{chrom}\t{a1}\t{a1+RES}\t{chrom}\t{a2}\t{a2+RES}\n")

for cls, rows in out.items():
    with open(f"ctcf_loops_{cls}_raw2fold.bedpe", "w") as fh:
        fh.writelines(rows)
    print(f"{cls:7s}: {len(rows):,}")
