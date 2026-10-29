"""
Transcription-quintile vs loop-change boxplots, varying the window around the
CTCF peak used to measure TT-seq log2FC: 500bp, 1kb, 2kb, 5kb, 10kb, 25kb (total
widths, centered on each peak).

Transcription: log2_int11pvsint1m.bw, mean over CTCF peaks (center +/- W/2) inside
               the 25 kb anchor bin; per loop = max of the two anchors.
Loop change  : log2(+aux/ctl) from SIP APA strengths in CTCFLoops.txt.
Loops        : strength<=20 set.
"""
import bisect, numpy as np, pyBigWig
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

BW="data/ttseq_bigwigs/log2_int11pvsint1m.bw"
BIN=25000; N_Q=5
WINDOWS=[("500bp",250),("1kb",500),("2kb",1000),("5kb",2500),("10kb",5000),("25kb",12500)]

def plain(c): return c[3:] if c.startswith("chr") else c
def resolve(bw,c):
    if bw.chroms(c) is not None: return c
    a=("chr"+c) if not c.startswith("chr") else c[3:]
    return a if bw.chroms(a) is not None else None

peaks={}
for l in open("data/loops/GSM3715311_CTCF_CT_peaks.bed"):
    if not l.strip(): continue
    c=l.split("\t"); peaks.setdefault(plain(c[0]),[]).append((int(c[1])+int(c[2]))//2)
for k in peaks: peaks[k].sort()

strength={}
with open("data/loops/ctl_totloops_strength.tab") as fh:
    next(fh)
    for c in (l.rstrip("\n").split("\t") for l in fh): strength[(c[1],int(c[2]),int(c[5]))]=float(c[7])
loops=[]
for l in open("data/loops/CTCFLoops.txt"):
    if not l.strip(): continue
    c=l.rstrip("\n").split("\t"); k=(c[0],int(c[1]),int(c[2]))
    if strength.get(k,0.0)<=20: loops.append((c[0],int(c[1]),int(c[2]),float(c[3]),float(c[5])))
loop_lfc=np.array([np.log2(pa/cl) for *_,cl,pa in loops])

bw=pyBigWig.open(BW)
# precompute, per anchor, the list of peak centers in its bin (shared across windows)
def peak_centers(ch,start):
    arr=peaks.get(plain(ch),[]); lo=bisect.bisect_left(arr,start); hi=bisect.bisect_left(arr,start+BIN)
    return arr[lo:hi]
anchor_pk={}
for ch,a1,a2,_,_ in loops:
    for a in (a1,a2):
        if (ch,a) not in anchor_pk: anchor_pk[(ch,a)]=peak_centers(ch,a)

def anchor_txn(ch,a,half):
    bwc=resolve(bw,plain(ch));
    if bwc is None: return np.nan
    clen=bw.chroms(bwc); vals=[]
    for mid in anchor_pk[(ch,a)]:
        s=max(0,mid-half); e=min(clen,mid+half)
        if s<e:
            st=bw.stats(bwc,s,e,type="mean")
            if st and st[0] is not None: vals.append(float(st[0]))
    return float(np.mean(vals)) if vals else np.nan

fig,axes=plt.subplots(2,3,figsize=(16,9)); axes=axes.ravel()
colors=plt.cm.coolwarm(np.linspace(0,1,N_Q))
print(f"loops (strength<=20) = {len(loops)}")
for w,(lab,half) in enumerate(WINDOWS):
    txn=np.array([np.nanmax([anchor_txn(ch,a1,half),anchor_txn(ch,a2,half)]) for ch,a1,a2,_,_ in loops])
    ok=np.isfinite(txn)&np.isfinite(loop_lfc)
    idx=np.where(ok)[0]; order=idx[np.argsort(txn[idx])]
    groups=np.array_split(order,N_Q)
    box=[loop_lfc[g] for g in groups]
    ax=axes[w]
    bp=ax.boxplot(box,showfliers=False,patch_artist=True,medianprops=dict(color="black"))
    for p,c in zip(bp["boxes"],colors): p.set_facecolor(c); p.set_alpha(0.85)
    ax.axhline(0,color="grey",lw=0.8,ls="--")
    ax.set_xticklabels([f"Q{i+1}\n({len(groups[i])})" for i in range(N_Q)],fontsize=7)
    ax.set_ylabel("loop change log2(+aux/ctl)",fontsize=8)
    u,p=mannwhitneyu(box[-1],box[0],alternative="less")
    meds=[np.median(b) for b in box]
    ax.set_title(f"peak ±{half}bp ({lab} window)\nQ1 med={meds[0]:+.2f} → Q5 med={meds[-1]:+.2f}  |  Q5<Q1 p={p:.1e}",
                 fontsize=8,fontweight="bold")
    ax.set_xlabel("TT-seq log2FC quintile (Q1 low → Q5 high)",fontsize=8)
    print(f"  {lab:>5} (±{half}): Q-medians={[round(m,3) for m in meds]}  Q5vsQ1 p={p:.2e}")
bw.close()
fig.suptitle("Loop change vs TT-seq transcription quintile — varying CTCF-peak window\n"
             "(transcription = log2_int11pvsint1m.bw at peaks; loop change = log2(+aux/ctl) SIP strength)",
             fontsize=12,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("quintile_windows_boxplots.png",dpi=120); fig.savefig("quintile_windows_boxplots.svg")
print("Saved quintile_windows_boxplots.png / .svg")
