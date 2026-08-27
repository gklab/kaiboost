#!/usr/bin/env python3
"""Decisive test: does the proband show more chromosome-level coverage
dispersion than a karyotypically normal control run through the same pipeline?

Variegated mosaic aneuploidy = several chromosomes each off by a few percent,
i.e. EXCESS chromosome-to-chromosome dispersion beyond technical noise. So the
statistic is the spread of per-chromosome coverage ratios, proband vs control.

Depth is equalised by binomial thinning of the deeper sample (repeated over
seeds), because dispersion shrinks with depth and the control is shallower.
Both samples go through the identical GC x repeat, leave-one-chromosome-out
model. chrX (1 copy, both males) is the truth anchor in each.

Control: GIAB HG002, PCR-free TruSeq, HiSeq 2500 2x148 — a different platform
from the proband's NovaSeq, so it bounds rather than exactly matches the
technical noise. That asymmetry is stated in the verdict.
"""
import glob
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
AUTO = [str(c) for c in range(1, 23)]
CHROMS = AUTO + ["X", "Y"]
N_SEEDS = 8


def load(pattern):
    counts = defaultdict(int)
    files = sorted(glob.glob(str(ROOT / pattern)))
    for fp in files:
        with open(fp) as f:
            f.readline()
            for line in f:
                c, b, n = line.split("\t")
                counts[(c, int(b))] += int(n)
    return counts, len(files)


def load_feats():
    feats = {}
    with open(ROOT / "cache" / "bin_features_100kb.tsv") as f:
        f.readline()
        for line in f:
            c, b, gc, rep = line.split("\t")
            feats[(c, int(b))] = (float(gc), float(rep))
    return feats


def design(gc, rep):
    g = (gc - 0.41) / 0.06
    r = (rep - 0.50) / 0.12
    return np.column_stack([np.ones_like(g), g, g**2, g**3, g**4,
                            r, r**2, r**3, g * r])


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


def chrom_ratios(obs, gc, rep, ch):
    """LOCO GC x repeat model -> median coverage ratio per chromosome."""
    med = np.median(obs[np.isin(ch, AUTO)])
    ok = (obs > 0.3 * med) & (obs < 2.0 * med)
    X, y = design(gc, rep), np.log(np.maximum(obs, 1))
    out = {}
    for c in CHROMS:
        train = ok & np.isin(ch, AUTO) & (ch != c)
        beta = robust_fit(X[train], y[train])
        sel = (ch == c) & ok
        if sel.sum() >= 20:
            out[c] = float(np.median(np.exp(y[sel] - X[sel] @ beta)))
    return out


def main():
    feats = load_feats()
    pro, n_pro = load("fastq_pass/L00*.bins.tsv")
    ctl, n_ctl = load("control_pass/ctrl_*.bins.tsv")
    tot_pro, tot_ctl = sum(pro.values()), sum(ctl.values())
    print(f"proband: {n_pro} units, {tot_pro/1e6:.1f}M binned reads")
    print(f"control: {n_ctl} units, {tot_ctl/1e6:.1f}M binned reads (HG002)")

    keys = sorted(set(pro) & set(ctl) & set(feats))
    keys = [k for k in keys if k[0] in CHROMS]
    gc = np.array([feats[k][0] for k in keys])
    rep = np.array([feats[k][1] for k in keys])
    ch = np.array([k[0] for k in keys])
    pro_c = np.array([pro[k] for k in keys], float)
    ctl_c = np.array([ctl[k] for k in keys], float)
    print(f"shared bins: {len(keys)}")

    p_thin = min(1.0, ctl_c.sum() / pro_c.sum())
    print(f"\nthinning proband to control depth (p={p_thin:.3f}), "
          f"{N_SEEDS} seeds\n")

    ctl_r = chrom_ratios(ctl_c, gc, rep, ch)
    pro_runs = []
    for s in range(N_SEEDS):
        rng = np.random.default_rng(1000 + s)
        thinned = rng.binomial(pro_c.astype(int), p_thin).astype(float)
        pro_runs.append(chrom_ratios(thinned, gc, rep, ch))
    pro_r = {c: float(np.mean([r[c] for r in pro_runs])) for c in ctl_r}

    print("chrom   proband(thinned)   control(HG002)    diff")
    for c in CHROMS:
        if c not in ctl_r or c not in pro_r:
            continue
        sd = np.std([r[c] for r in pro_runs])
        print(f"{c:>5}   {pro_r[c]:.4f} +-{sd:.4f}      {ctl_r[c]:.4f}      "
              f"{pro_r[c]-ctl_r[c]:+.4f}")

    pa = np.array([pro_r[c] for c in AUTO if c in pro_r])
    ca = np.array([ctl_r[c] for c in AUTO if c in ctl_r])
    sd_pro = float(pa.std())
    sd_ctl = float(ca.std())
    per_seed_sd = [np.std([r[c] for c in AUTO if c in r]) for r in pro_runs]

    print(f"\n=== autosomal dispersion (the aneuploidy statistic) ===")
    print(f"  proband  SD = {sd_pro*100:.2f}%  (per-seed range "
          f"{min(per_seed_sd)*100:.2f}-{max(per_seed_sd)*100:.2f}%)")
    print(f"  control  SD = {sd_ctl*100:.2f}%")
    print(f"  ratio    = {sd_pro/sd_ctl:.2f}x")

    print(f"\n=== truth anchor (chrX = 1 copy in both males) ===")
    for name, r in (("proband", pro_r), ("control", ctl_r)):
        if "X" in r:
            print(f"  {name}: chrX {r['X']:.4f} vs 0.5000 -> error "
                  f"{abs(r['X']/0.5-1)*100:.2f}%")

    if len(pa) == len(ca):
        cc = float(np.corrcoef(pa, ca)[0, 1])
        print(f"\n=== shared-artifact check ===")
        print(f"  corr(proband deviations, control deviations) = {cc:+.3f}")
        print("  (high corr => common pipeline artifact; ~0 => sample-specific"
              " noise, not correctable by subtracting a control)")

    print("\n=== VERDICT ===")
    if sd_pro <= sd_ctl * 1.2:
        print("  NEGATIVE. The proband's chromosome-level spread does not exceed")
        print("  a karyotypically normal sample's at matched depth. No mosaic")
        print("  aneuploidy is detectable by this method; the earlier"
              " chr16/17/18/19/21 'gains' are within technical noise.")
    else:
        print(f"  EXCESS DISPERSION ({sd_pro/sd_ctl:.2f}x control). Not proof —")
        print("  control is a different platform — but worth deeper follow-up.")
    print("\n  Limitation: control is HiSeq 2500 2x148 PCR-free vs proband")
    print("  NovaSeq; platform differences are NOT controlled. This bounds the")
    print("  noise floor, it does not certify the absence of low-level mosaicism.")


if __name__ == "__main__":
    main()
