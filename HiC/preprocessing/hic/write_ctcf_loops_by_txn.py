"""
Write CTCF loops split into transcription-change quintiles (bedpe), for SIPmeta.

Each loop gets one transcription value = MAX of its two anchors' peak-refined
TT-seq log2FC (CTCF peaks center +/-2500bp inside the 25 kb anchor bin).
Loops are split into 5 equal-count quintiles (Q1 lowest -> Q5 highest txn) and
written as:  chr  a1  a1+25kb  chr  a2  a2+25kb   (chr-prefixed, tab-separated).
"""
import bisect
import numpy as np
import pyBigWig

BW_PATH   = "data/ttseq_bigwigs/log2_int11pvsint1m.bw"
LOOP_PATH = "data/loops/CTCFLoops.txt"
PEAK_PATH = "data/loops/GSM3715311_CTCF_CT_peaks.bed"
BIN       = 25000
PEAK_WIN  = 2500
N_Q       = 5
OUT_STEM  = "ctcf_loops_txnQ"          # -> ctcf_loops_txnQ1_int11vsctl.bedpe ...
OUT_SUFF  = "_int11vsctl.bedpe"


def plain(c):
    return c[3:] if c.startswith("chr") else c


def resolve_chrom(bw, chrom):
    if bw.chroms(chrom) is not None:
        return chrom
    alt = ("chr" + chrom) if not chrom.startswith("chr") else chrom[3:]
    return alt if bw.chroms(alt) is not None else None


def load_peaks():
    d = {}
    for line in open(PEAK_PATH):
        if not line.strip():
            continue
        c = line.split("\t")
        d.setdefault(plain(c[0]), []).append((int(c[1]) + int(c[2])) // 2)
    for ch in d:
        d[ch].sort()
    return d


def main():
    peaks = load_peaks()
    bw = pyBigWig.open(BW_PATH)
    cache = {}

    def anchor_txn(chrom_plain, start):
        key = (chrom_plain, start)
        if key in cache:
            return cache[key]
        bwc = resolve_chrom(bw, chrom_plain)
        vals = []
        arr = peaks.get(chrom_plain, [])
        lo = bisect.bisect_left(arr, start); hi = bisect.bisect_left(arr, start + BIN)
        if bwc is not None:
            clen = bw.chroms(bwc)
            for mid in arr[lo:hi]:
                s = max(0, mid - PEAK_WIN); e = min(clen, mid + PEAK_WIN)
                if s < e:
                    st = bw.stats(bwc, s, e, type="mean")
                    if st and st[0] is not None:
                        vals.append(float(st[0]))
        val = float(np.mean(vals)) if vals else np.nan
        cache[key] = val
        return val

    loops = []   # (chrom_orig, a1, a2, txn)
    for line in open(LOOP_PATH):
        if not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        chrom = c[0]; a1, a2 = int(c[1]), int(c[2])
        t1 = anchor_txn(plain(chrom), a1)
        t2 = anchor_txn(plain(chrom), a2)
        txn = np.nanmax([t1, t2]) if (np.isfinite(t1) or np.isfinite(t2)) else np.nan
        loops.append((chrom, a1, a2, txn))
    bw.close()

    valid = [l for l in loops if np.isfinite(l[3])]
    dropped = len(loops) - len(valid)
    valid.sort(key=lambda l: l[3])
    groups = np.array_split(np.arange(len(valid)), N_Q)

    print(f"CTCF loops {len(loops):,} | finite txn {len(valid):,} | dropped (no peak/data) {dropped}")
    for q, gi in enumerate(groups):
        members = [valid[k] for k in gi]
        fn = f"{OUT_STEM}{q+1}{OUT_SUFF}"
        with open(fn, "w") as fh:
            for chrom, a1, a2, txn in members:
                fh.write(f"{chrom}\t{a1}\t{a1 + BIN}\t{chrom}\t{a2}\t{a2 + BIN}\n")
        tlo, thi = members[0][3], members[-1][3]
        print(f"  Q{q+1}: {len(members):,} loops  txn log2FC [{tlo:.3f}, {thi:.3f}]  -> {fn}")


if __name__ == "__main__":
    main()
