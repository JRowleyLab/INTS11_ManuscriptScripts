"""
Relate the degree (magnitude) and shape (spatial profile) of INTS11-induced TT-seq
transcription changes to the corresponding CRUSH (Hi-C compartment) changes around TSSs.

Uses the matched, strand-aware, TSS-anchored matrices (same 8,874 genes, same order,
20 kb up / 4 kb down, 50 bp bins -> 480 bins, TSS at bin 400) for two pairings:

  control : ttseq_int11pvsint1m_stranded.npz   vs  hic_delta_stranded.npz        (INT11vsctl)
  noaux   : ttseq_int11pvsint11m_stranded.npz  vs  hic_delta_noaux_stranded.npz  (INT11vsnoaux)

Outputs three figure types per pairing + a stats summary. No matrix recomputation.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib import cm
from scipy.stats import spearmanr, pearsonr

# ── geometry (mirrors plot_ttseq_delta_heatmap.py) ─────────────────────────────
UPSTREAM, DOWNSTREAM, BINSIZE = 20000, 4000, 50
N_BINS  = (UPSTREAM + DOWNSTREAM) // BINSIZE   # 480
TSS_BIN = UPSTREAM // BINSIZE                  # 400
PROX_LO = TSS_BIN - 40                         # -2 kb
PROX_HI = N_BINS                               # +4 kb

X_TICKS  = [0, 100, 200, 300, TSS_BIN, TSS_BIN + 40, TSS_BIN + 80]
X_LABELS = ["-20kb", "-15kb", "-10kb", "-5kb", "TSS", "+2kb", "+4kb"]

PAIRINGS = {
    "control": dict(
        tt="ttseq_int11pvsint1m_stranded.npz",  crush="hic_delta_stranded.npz",
        tt_lab="TT-seq log2FC (INTS11 vs INTS1)",
        crush_lab="CRUSH Δ (INT11 vs ctl)"),
    "noaux": dict(
        tt="ttseq_int11pvsint11m_stranded.npz", crush="hic_delta_noaux_stranded.npz",
        tt_lab="TT-seq log2FC (INTS11+ vs INTS11-)",
        crush_lab="CRUSH Δ (INT11 vs noAux)"),
}

div_cmap = mcolors.LinearSegmentedColormap.from_list(
    "bwo", ["#4393C3", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#E08214"], N=512)


# ── helpers ────────────────────────────────────────────────────────────────────
def load_pair(cfg):
    tt = np.load(cfg["tt"], allow_pickle=True)["matrix"].astype(np.float64)
    cr = np.load(cfg["crush"], allow_pickle=True)["matrix"].astype(np.float64)
    return tt, cr


def summarize(mat, lo, hi):
    """Per-gene mean over bin window [lo:hi], ignoring NaN."""
    return np.nanmean(mat[:, lo:hi], axis=1)


def winsor(a, lo=1, hi=99):
    ql, qh = np.nanpercentile(a, [lo, hi])
    return np.clip(a, ql, qh)


def finite_pair(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m], int(m.sum())


# ── Analysis 1: per-gene magnitude scatter + correlation ───────────────────────
def scatter_panel(ax, s_tt, s_cr, cfg):
    x, y, n = finite_pair(s_tt, s_cr)
    rho, p_s = spearmanr(x, y)
    r,   p_p = pearsonr(x, y)
    xl = np.nanpercentile(x, [1, 99]); yl = np.nanpercentile(y, [1, 99])
    ax.hexbin(x, y, gridsize=60, bins="log", cmap="viridis",
              extent=[xl[0], xl[1], yl[0], yl[1]], mincnt=1)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlim(*xl); ax.set_ylim(*yl)
    ax.set_xlabel(f"Per-gene {cfg['tt_lab']}\n(mean -2kb→+4kb)", fontsize=8)
    ax.set_ylabel(f"Per-gene {cfg['crush_lab']}\n(mean -2kb→+4kb)", fontsize=8)
    ax.set_title(f"Spearman ρ = {rho:.3f}  (p={p_s:.1e})\n"
                 f"Pearson r = {r:.3f}  (p={p_p:.1e})  |  n = {n:,}",
                 fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=7)
    return dict(rho=rho, p_s=p_s, r=r, p_p=p_p, n=n)


def make_scatter(per_gene):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    stats = {}
    for ax, tag in zip(axes, PAIRINGS):
        s_tt, s_cr = per_gene[tag]
        stats[tag] = scatter_panel(ax, s_tt, winsor(s_cr), PAIRINGS[tag])
        ax.text(0.02, 0.97, tag, transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="top", ha="left",
                bbox=dict(fc="white", ec="grey", alpha=0.7, pad=2))
    fig.suptitle("Degree of TT-seq change vs CRUSH change (per gene, TSS-proximal)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("ttseq_crush_scatter.png", dpi=150, bbox_inches="tight")
    fig.savefig("ttseq_crush_scatter.svg", bbox_inches="tight")
    plt.close(fig)
    print("Saved ttseq_crush_scatter.png / .svg")
    return stats


# ── Analysis 2: quantile metagene profiles (signed + magnitude) ────────────────
def winsor_mat(mat, lo=1, hi=99):
    """Clip a matrix to global percentiles (for robust CRUSH metagene means)."""
    ql, qh = np.nanpercentile(mat, [lo, hi])
    return np.clip(mat, ql, qh)


def quantile_profiles(ax, rank_by, profile_mat, n_q, title, ylab, mode):
    """Split genes into n_q quantiles by rank_by; plot mean profile per group.

    mode='signed'   : rank by signed value, coolwarm Q1(most ↓)…Q5(most ↑)
    mode='magnitude': rank by |value|,      viridis  Q1(smallest |Δ|)…Q5(largest |Δ|)
    """
    if mode == "signed":
        key = rank_by; cmap = cm.coolwarm; lo_lab, hi_lab = "most ↓", "most ↑"
    else:
        key = np.abs(rank_by); cmap = cm.viridis; lo_lab, hi_lab = "smallest |Δ|", "largest |Δ|"
    valid = np.isfinite(key)
    idx = np.where(valid)[0][np.argsort(key[valid])]
    groups = np.array_split(idx, n_q)
    colors = cmap(np.linspace(0, 1, n_q))
    x = np.arange(N_BINS)
    for q, (g, c) in enumerate(zip(groups, colors)):
        prof = np.nanmean(profile_mat[g], axis=0)
        lab = f"Q{q+1}"
        if q == 0:        lab += f" ({lo_lab})"
        elif q == n_q-1:  lab += f" ({hi_lab})"
        ax.plot(x, prof, color=c, lw=1.3, label=lab)
    ax.axvline(TSS_BIN, color="black", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlim(0, N_BINS)
    ax.set_xticks(X_TICKS); ax.set_xticklabels(X_LABELS, fontsize=7)
    ax.set_ylabel(ylab, fontsize=8)
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6.5, loc="upper left", frameon=False, ncol=1)


def make_quantile(tag, tt, cr, s_tt, s_cr, cfg, n_q=5):
    """2x2: rows = signed / magnitude ranking; cols = CRUSH-by-TTseq, TTseq-by-CRUSH."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    cr_w = winsor_mat(cr)
    for row, mode in enumerate(("signed", "magnitude")):
        quantile_profiles(
            axes[row, 0], s_tt, cr_w, n_q,
            f"CRUSH metagene · genes ranked by {mode} {cfg['tt_lab']}",
            cfg["crush_lab"], mode)
        quantile_profiles(
            axes[row, 1], s_cr, tt, n_q,
            f"TT-seq metagene · genes ranked by {mode} {cfg['crush_lab']}",
            cfg["tt_lab"], mode)
    fig.suptitle(f"Quantile metagenes (TSS-proximal -2kb→+4kb ranking)  [{tag}]\n"
                 "top row = signed ranking,  bottom row = |change| (magnitude) ranking",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"ttseq_crush_quantile_{tag}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"ttseq_crush_quantile_{tag}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved ttseq_crush_quantile_{tag}.png / .svg")


# ── Analysis 2c: one shared ranking (signed TT-seq) for both signals ───────────
def make_signed_shared(tag, tt, cr, s_tt, cfg, n_q=5):
    """Rank genes ONCE by signed TT-seq log2FC; plot TT-seq and CRUSH metagenes
    on that single shared ranking (same gene groups in both panels)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    # both panels ranked by the SAME signal (s_tt) -> identical quintile groups
    quantile_profiles(
        axes[0], s_tt, tt, n_q,
        f"TT-seq metagene  (ranking signal)",
        cfg["tt_lab"], "signed")
    quantile_profiles(
        axes[1], s_tt, winsor_mat(cr), n_q,
        f"CRUSH metagene  (same TT-seq ranking)",
        cfg["crush_lab"], "signed")
    fig.suptitle(f"Signed TT-seq ranking — both signals on one ordering  [{tag}]\n"
                 f"genes ranked by signed {cfg['tt_lab']} (mean -2kb→+4kb), "
                 "Q1 = most ↓ … Q5 = most ↑",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(f"ttseq_crush_signedrank_{tag}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"ttseq_crush_signedrank_{tag}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved ttseq_crush_signedrank_{tag}.png / .svg")


# ── Analysis 2b: up / down / unchanged class metagenes ─────────────────────────
def classify(s):
    """Sign classes with a neutral zone gated at the median |change|."""
    gate = np.nanmedian(np.abs(s))
    cls = np.full(s.shape, "none", dtype=object)
    cls[np.isfinite(s) & (s >  gate)] = "up"
    cls[np.isfinite(s) & (s < -gate)] = "down"
    cls[np.isfinite(s) & (np.abs(s) <= gate)] = "unchanged"
    return cls, gate


def class_panel(ax, classifier, profile_mat, title, ylab):
    cls, gate = classify(classifier)
    style = dict(up=("firebrick", "up"), down=("royalblue", "down"),
                 unchanged=("grey", "unchanged"))
    x = np.arange(N_BINS)
    counts = {}
    for key, (color, name) in style.items():
        g = np.where(cls == key)[0]
        counts[key] = len(g)
        prof = np.nanmean(profile_mat[g], axis=0)
        ax.plot(x, prof, color=color, lw=1.5, label=f"{name} (n={len(g):,})")
    ax.axvline(TSS_BIN, color="black", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlim(0, N_BINS)
    ax.set_xticks(X_TICKS); ax.set_xticklabels(X_LABELS, fontsize=7)
    ax.set_ylabel(ylab, fontsize=8)
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    return counts


def make_classes(tag, tt, cr, s_tt, s_cr, cfg):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    c1 = class_panel(axes[0], s_tt, winsor_mat(cr),
                     f"CRUSH metagene by TT-seq class\n(up/down/unchanged by {cfg['tt_lab']})",
                     cfg["crush_lab"])
    c2 = class_panel(axes[1], s_cr, tt,
                     f"TT-seq metagene by CRUSH class\n(up/down/unchanged by {cfg['crush_lab']})",
                     cfg["tt_lab"])
    fig.suptitle(f"Up / down / unchanged class metagenes  [{tag}]\n"
                 "classes: sign of TSS-proximal change, neutral zone = |Δ| ≤ median(|Δ|)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(f"ttseq_crush_classes_{tag}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"ttseq_crush_classes_{tag}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved ttseq_crush_classes_{tag}.png / .svg")
    return dict(by_ttseq=c1, by_crush=c2)


# ── Analysis 3: bin-level 2D density + spatial cross-correlation ───────────────
def make_binlevel(tag, tt, cr, cfg):
    fig = plt.figure(figsize=(14, 4.5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1.1, 1], wspace=0.32)

    # (a) 2D hexbin of all matched (gene,bin) pairs
    ax_a = fig.add_subplot(gs[0])
    xf, yf, n = finite_pair(tt.ravel(), cr.ravel())
    rho, _ = spearmanr(xf, yf)
    r,   _ = pearsonr(xf, yf)
    xl = np.nanpercentile(xf, [1, 99]); yl = np.nanpercentile(yf, [1, 99])
    ax_a.hexbin(xf, yf, gridsize=80, bins="log", cmap="magma",
                extent=[xl[0], xl[1], yl[0], yl[1]], mincnt=1)
    ax_a.axhline(0, color="white", lw=0.4, ls="--")
    ax_a.axvline(0, color="white", lw=0.4, ls="--")
    ax_a.set_xlim(*xl); ax_a.set_ylim(*yl)
    ax_a.set_xlabel(cfg["tt_lab"], fontsize=8)
    ax_a.set_ylabel(cfg["crush_lab"], fontsize=8)
    ax_a.set_title(f"(a) Bin-level density\nρ={rho:.3f}  r={r:.3f}  (n={n:,})",
                   fontsize=8, fontweight="bold")
    ax_a.tick_params(labelsize=7)

    # (b) per-position Spearman across genes
    ax_b = fig.add_subplot(gs[1])
    pos_rho = np.full(N_BINS, np.nan)
    for b in range(N_BINS):
        x = tt[:, b]; y = cr[:, b]
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() > 10:
            pos_rho[b] = spearmanr(x[m], y[m]).correlation
    xaxis = np.arange(N_BINS)
    ax_b.fill_between(xaxis, 0, pos_rho, color="purple", alpha=0.35)
    ax_b.plot(xaxis, pos_rho, color="purple", lw=1.1)
    ax_b.axvline(TSS_BIN, color="black", lw=0.8)
    ax_b.axhline(0, color="grey", lw=0.5, ls="--")
    ax_b.set_xlim(0, N_BINS); ax_b.set_xticks(X_TICKS)
    ax_b.set_xticklabels(X_LABELS, fontsize=7)
    ax_b.set_ylabel("Spearman ρ (TT-seq vs CRUSH)\nacross genes, per position", fontsize=8)
    ax_b.set_title("(b) Per-position correlation", fontsize=8, fontweight="bold")
    ax_b.tick_params(labelsize=7)

    # (c) spatial cross-correlation of the z-scored mean metagene profiles
    ax_c = fig.add_subplot(gs[2])
    p_tt = np.nanmean(tt, axis=0); p_cr = np.nanmean(cr, axis=0)
    zt = (p_tt - p_tt.mean()) / p_tt.std()
    zc = (p_cr - p_cr.mean()) / p_cr.std()
    full = np.correlate(zt, zc, mode="full") / len(zt)
    lags = np.arange(-(N_BINS - 1), N_BINS)
    win = np.abs(lags) <= 100
    lags_w, cc_w = lags[win], full[win]
    peak = lags_w[np.argmax(cc_w)]
    ax_c.axvline(0, color="grey", lw=0.5, ls="--")
    ax_c.plot(lags_w * BINSIZE / 1000.0, cc_w, color="teal", lw=1.2)
    ax_c.axvline(peak * BINSIZE / 1000.0, color="red", lw=0.8, ls=":")
    ax_c.set_xlabel("Lag (kb)  [TT-seq shifted vs CRUSH]", fontsize=8)
    ax_c.set_ylabel("Cross-correlation\n(z-scored mean profiles)", fontsize=8)
    ax_c.set_title(f"(c) Spatial cross-correlation\npeak lag = {peak*BINSIZE/1000:.2f} kb",
                   fontsize=8, fontweight="bold")
    ax_c.tick_params(labelsize=7)

    fig.suptitle(f"Shape relationship: TT-seq vs CRUSH  [{tag}]",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"ttseq_crush_binlevel_{tag}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"ttseq_crush_binlevel_{tag}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved ttseq_crush_binlevel_{tag}.png / .svg")
    return dict(bin_rho=rho, bin_r=r, bin_n=n, peak_lag_bins=int(peak),
                peak_lag_kb=peak * BINSIZE / 1000.0)


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    per_gene = {}
    summ = {}
    for tag, cfg in PAIRINGS.items():
        print(f"\n=== {tag} ===")
        tt, cr = load_pair(cfg)
        s_tt_prox = summarize(tt, PROX_LO, PROX_HI)
        s_cr_prox = summarize(cr, PROX_LO, PROX_HI)
        s_tt_full = summarize(tt, 0, N_BINS)
        s_cr_full = summarize(cr, 0, N_BINS)
        per_gene[tag] = (s_tt_prox, s_cr_prox)

        # window correlations (report both)
        rec = {}
        for wname, (a, b) in dict(proximal=(s_tt_prox, s_cr_prox),
                                  full=(s_tt_full, s_cr_full)).items():
            x, y, n = finite_pair(a, winsor(b))
            rec[wname] = dict(rho=spearmanr(x, y).correlation,
                              r=pearsonr(x, y)[0], n=n)
        summ[tag] = rec

        make_quantile(tag, tt, cr, s_tt_prox, s_cr_prox, cfg)
        make_signed_shared(tag, tt, cr, s_tt_prox, cfg)
        summ[tag]["classes"] = make_classes(tag, tt, cr, s_tt_prox, s_cr_prox, cfg)
        summ[tag]["binlevel"] = make_binlevel(tag, tt, cr, cfg)

    scat = make_scatter(per_gene)
    for tag in PAIRINGS:
        summ[tag]["scatter_proximal"] = scat[tag]

    # ── stats summary file ──────────────────────────────────────────────────────
    with open("ttseq_crush_stats.txt", "w") as fh:
        fh.write("TT-seq vs CRUSH relationship — summary statistics\n")
        fh.write("=" * 60 + "\n")
        for tag, cfg in PAIRINGS.items():
            r = summ[tag]
            fh.write(f"\n[{tag}]  {cfg['tt_lab']}  vs  {cfg['crush_lab']}\n")
            fh.write(f"  Per-gene correlation (TSS-proximal -2kb->+4kb): "
                     f"Spearman rho={r['proximal']['rho']:.3f}  "
                     f"Pearson r={r['proximal']['r']:.3f}  n={r['proximal']['n']:,}\n")
            fh.write(f"  Per-gene correlation (full -20kb->+4kb window):  "
                     f"Spearman rho={r['full']['rho']:.3f}  "
                     f"Pearson r={r['full']['r']:.3f}  n={r['full']['n']:,}\n")
            bl = r["binlevel"]
            fh.write(f"  Bin-level (all gene x bin) correlation:          "
                     f"Spearman rho={bl['bin_rho']:.3f}  "
                     f"Pearson r={bl['bin_r']:.3f}  n={bl['bin_n']:,}\n")
            fh.write(f"  Spatial cross-correlation peak lag:              "
                     f"{bl['peak_lag_kb']:.2f} kb ({bl['peak_lag_bins']} bins)\n")
            cl = r["classes"]
            fh.write(f"  Class counts by TT-seq  (up/unch/down): "
                     f"{cl['by_ttseq']['up']:,} / {cl['by_ttseq']['unchanged']:,} / "
                     f"{cl['by_ttseq']['down']:,}\n")
            fh.write(f"  Class counts by CRUSH   (up/unch/down): "
                     f"{cl['by_crush']['up']:,} / {cl['by_crush']['unchanged']:,} / "
                     f"{cl['by_crush']['down']:,}\n")
    print("\nSaved ttseq_crush_stats.txt")
    print(open("ttseq_crush_stats.txt").read())


if __name__ == "__main__":
    main()
