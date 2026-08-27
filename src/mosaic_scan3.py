#!/usr/bin/env python3
"""Mappability-aware mosaic aneuploidy scan (supersedes mosaic_scan2.py).

mosaic_scan2 corrected coverage for GC only and flagged +1.5-4% "gains" on
chr16/17/18/19/21/22 — but chrX (known 1 copy in this male) came out 1.6%
off its true 0.5, i.e. the GC-only model has a demonstrated systematic error
of the same order as the flagged signal, and every flagged 10 Mb segment sat
in pericentromeric / segmental-duplication sequence. Both point at
between-chromosome mappability differences, not copy number.

This version therefore models coverage as a smooth function of BOTH GC and
repeat content (lowercase fraction of the soft-masked UCSC reference, a free
RepeatMasker proxy for mappability), fitted robustly in log space.

Self-absorption guard: for each chromosome the model is fitted LEAVE-ONE-
CHROMOSOME-OUT, so a chromosome sitting in its own corner of feature space
(chr19 = high GC + Alu-rich) cannot explain away its own deviation.

Truth anchors: chrX and chrY (1 copy each, male) must land at ratio 0.5 —
their residual error IS the method's empirical detection floor, and any
autosomal claim must clear it.
"""
import glob
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

ROOT = Path(__file__).resolve().parent.parent
BIN = 100_000
AUTO = [str(c) for c in range(1, 23)]
CHROMS = AUTO + ["X", "Y"]
FEAT = ROOT / "cache" / "bin_features_100kb.tsv"
rng = np.random.default_rng(20260827)


PATTERN = sys.argv[1] if len(sys.argv) > 1 else "fastq_pass/L00*.bins.tsv"
LABEL = sys.argv[2] if len(sys.argv) > 2 else "mosaic_scan3"


def load_counts():
    counts = defaultdict(int)
    files = sorted(glob.glob(str(ROOT / PATTERN)))
    if not files:
        sys.exit("no bins files found")
    for fp in files:
        with open(fp) as f:
            f.readline()
            for line in f:
                c, b, n = line.split("\t")
                counts[(c, int(b))] += int(n)
    print(f"[bins] {len(files)} lanes, {len(counts)} bins, "
          f"{sum(counts.values())/1e6:.1f}M reads")
    return counts


def load_features(needed):
    """(chrom, bin) -> (gc, repeat_frac). Cached; one FASTA pass on first run."""
    if FEAT.exists():
        feats = {}
        with open(FEAT) as f:
            f.readline()
            for line in f:
                c, b, gc, rep = line.split("\t")
                feats[(c, int(b))] = (float(gc), float(rep))
        print(f"[feat] loaded {len(feats)} cached bins")
        return feats

    ref = sorted(glob.glob(str(ROOT / "cache" / "*.fa")),
                 key=lambda p: -Path(p).stat().st_size)[0]
    print(f"[feat] computing from {ref} (one-time)")
    fa = Fasta(ref, rebuild=False)
    key = {c.removeprefix("chr"): c for c in fa.keys()}
    by_chrom = defaultdict(list)
    for c, b in needed:
        by_chrom[c].append(b)
    feats = {}
    n_lower_total = 0
    for c in CHROMS:
        if c not in key:
            continue
        c_fa, seq_len = key[c], len(fa[key[c]])
        for b in sorted(by_chrom[c]):
            s = fa[c_fa][b * BIN:min((b + 1) * BIN, seq_len)].seq  # keep case
            if not s:
                continue
            n_n = s.count("N") + s.count("n")
            if n_n > 0.1 * len(s):
                continue
            acgt = len(s) - n_n
            if acgt < 0.5 * BIN:
                continue
            lower = sum(1 for ch in s if ch.islower())
            n_lower_total += lower
            up = s.upper()
            gc = (up.count("G") + up.count("C")) / acgt
            feats[(c, b)] = (gc, lower / acgt)
        print(f"  chr{c} done", flush=True)
    if n_lower_total == 0:
        sys.exit("reference is NOT soft-masked — repeat proxy unavailable")
    with open(FEAT, "w") as out:
        out.write("chrom\tbin\tgc\trepeat_frac\n")
        for (c, b), (gc, rep) in sorted(feats.items()):
            out.write(f"{c}\t{b}\t{gc:.5f}\t{rep:.5f}\n")
    print(f"[feat] wrote {len(feats)} bins -> {FEAT}")
    return feats


def design(gc, rep):
    """Polynomial basis in standardised GC and repeat fraction, plus interaction."""
    g = (gc - 0.41) / 0.06
    r = (rep - 0.50) / 0.12
    cols = [np.ones_like(g), g, g**2, g**3, g**4, r, r**2, r**3, g * r]
    return np.column_stack(cols)


def robust_fit(X, y, iters=4):
    w = np.ones(len(y), bool)
    beta = None
    for _ in range(iters):
        beta, *_ = np.linalg.lstsq(X[w], y[w], rcond=None)
        res = y - X @ beta
        c = np.median(res[w])
        mad = np.median(np.abs(res[w] - c)) * 1.4826 or 1e-9
        w = np.abs(res - c) < 3 * mad
    return beta


def main():
    counts = load_counts()
    feats = load_features(set(counts))

    keys = [k for k in counts if k in feats and k[0] in CHROMS]
    obs = np.array([counts[k] for k in keys], float)
    gc = np.array([feats[k][0] for k in keys])
    rep = np.array([feats[k][1] for k in keys])
    ch = np.array([k[0] for k in keys])
    bn = np.array([k[1] for k in keys])

    med = np.median(obs[np.isin(ch, AUTO)])
    usable = (obs > 0.3 * med) & (obs < 2.0 * med)
    print(f"[fit] {usable.sum()}/{len(obs)} bins usable "
          f"(gc {gc.min():.2f}-{gc.max():.2f}, repeat {rep.min():.2f}-{rep.max():.2f})")

    X = design(gc, rep)
    y = np.log(np.maximum(obs, 1))
    ratio = np.full(len(obs), np.nan)

    for c in CHROMS:                       # leave-one-chromosome-out fitting
        train = usable & np.isin(ch, AUTO) & (ch != c)
        beta = robust_fit(X[train], y[train])
        sel = ch == c
        ratio[sel] = np.exp(y[sel] - X[sel] @ beta)

    print("\nchrom  n_bins   ratio    ci95         implied_f   mean_gc mean_rep  verdict")
    rows = []
    for c in CHROMS:
        sel = (ch == c) & usable & ~np.isnan(ratio)
        r = ratio[sel]
        if len(r) < 20:
            continue
        m = float(np.median(r))
        boots = np.array([np.median(rng.choice(r, len(r))) for _ in range(400)])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        exp0 = 0.5 if c in ("X", "Y") else 1.0
        f = 2 * (m - exp0) / exp0
        sig = "" if lo <= exp0 <= hi else ("GAIN?" if m > exp0 else "LOSS?")
        rows.append((c, len(r), m, lo, hi, f, float(gc[sel].mean()),
                     float(rep[sel].mean()), sig))
        print(f"{c:>5} {len(r):7d}  {m:.4f}  {lo:.4f}-{hi:.4f}  {f:+8.3f}   "
              f"{gc[sel].mean():.3f}   {rep[sel].mean():.3f}   {sig}")

    with open(ROOT / "results" / f"{LABEL}.tsv", "w") as out:
        out.write("chrom\tn_bins\tratio\tci_lo\tci_hi\timplied_f\tmean_gc\t"
                  "mean_repeat\tverdict\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # --- empirical detection floor from the two known-truth chromosomes ---
    print("\n[truth anchors] X and Y are 1 copy in this male; deviation from"
          " 0.5 = method error")
    floor = 0.0
    for c in ("X", "Y"):
        row = next((r for r in rows if r[0] == c), None)
        if row:
            err = abs(row[2] / 0.5 - 1)
            floor = max(floor, err)
            print(f"  chr{c}: ratio {row[2]:.4f} vs 0.5000 -> relative error "
                  f"{err*100:.2f}% (= {err*200:.1f}% apparent mosaic fraction)")
    print(f"  => empirical floor: |implied_f| must exceed ~{floor*200:.1f}% to mean anything")

    print("\n[verdict] autosomes clearing the empirical floor:")
    hits = [r for r in rows if r[0] in AUTO and abs(r[5]) > floor * 2]
    if not hits:
        print("  NONE — no autosomal mosaic aneuploidy is detectable above the"
              " method's own systematic error")
    for r in hits:
        print(f"  chr{r[0]}: implied_f {r[5]:+.3f} (ratio {r[2]:.4f}, "
              f"CI {r[3]:.4f}-{r[4]:.4f})")

    # does any residual chromosome effect still track repeat content?
    au = [r for r in rows if r[0] in AUTO]
    rr = np.array([r[2] for r in au])
    pr = np.array([r[7] for r in au])
    pg = np.array([r[6] for r in au])
    print(f"\n[residual bias] corr(chrom ratio, mean repeat) = "
          f"{np.corrcoef(rr, pr)[0,1]:+.3f}; corr(ratio, mean GC) = "
          f"{np.corrcoef(rr, pg)[0,1]:+.3f}  (|r|>0.5 => residual artifact)")
    print(f"\nwrote results/{LABEL}.tsv")


if __name__ == "__main__":
    main()
