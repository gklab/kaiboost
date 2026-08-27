#!/usr/bin/env python3
"""GC-corrected mosaic aneuploidy scan from the FASTQ sample-pass bin counts.

Inputs: fastq_pass/L00*.bins.tsv (summed) + cache/GRCh38_primary.fa (GC per bin).
Model: expected count per bin = median count of its 1%-GC stratum (fit on
autosomes only); ratio = observed/expected. Per-chromosome robust mean ratio
with a bin-level bootstrap CI; implied mosaic fraction f = 2*|ratio-1|.
Also flags 10 Mb segments deviating from their own chromosome (segmental CNV).
Sex chromosomes reported against expected ratio 0.5 (male).
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
rng = np.random.default_rng(20260827)


def load_bins():
    counts = defaultdict(int)
    files = sorted(glob.glob(str(ROOT / "fastq_pass" / "L00*.bins.tsv")))
    if not files:
        sys.exit("no bins files found")
    print(f"[bins] merging {len(files)} lane files")
    for fp in files:
        with open(fp) as f:
            f.readline()
            for line in f:
                c, b, n = line.split("\t")
                counts[(c, int(b))] += int(n)
    return counts


def gc_per_bin(needed):
    import glob as g
    ref = sorted(g.glob(str(ROOT / "cache" / "*.fa")), key=lambda p: -Path(p).stat().st_size)[0]
    print(f"[gc] reference: {ref}")
    fa = Fasta(ref, rebuild=False)
    key = {c.removeprefix("chr"): c for c in fa.keys()}  # tolerate chr-prefixed refs
    gc = {}
    by_chrom = defaultdict(list)
    for c, b in needed:
        by_chrom[c].append(b)
    for c, bs in by_chrom.items():
        if c not in key:
            continue
        c_fa = key[c]
        seq_len = len(fa[c_fa])
        for b in bs:
            s = str(fa[c_fa][b * BIN:min((b + 1) * BIN, seq_len)]).upper()
            if not s:
                continue
            n_n = s.count("N")
            if n_n > 0.1 * len(s):
                continue
            acgt = len(s) - n_n
            gc[(c, b)] = (s.count("G") + s.count("C")) / max(acgt, 1)
        print(f"  gc: chr{c} done", flush=True)
    return gc


def main():
    counts = load_bins()
    print(f"[bins] {len(counts)} bins, total reads {sum(counts.values())/1e6:.1f}M")
    gc = gc_per_bin(set(counts))

    keys = [k for k in counts if k in gc]
    obs = np.array([counts[k] for k in keys], dtype=float)
    gcs = np.array([gc[k] for k in keys])
    chroms = np.array([k[0] for k in keys])
    bins_ = np.array([k[1] for k in keys])

    # fit expected-vs-GC on autosomes, trimming coverage outliers
    is_auto = np.isin(chroms, AUTO)
    med = np.median(obs[is_auto])
    ok = is_auto & (obs > 0.3 * med) & (obs < 2.0 * med)
    strata = np.clip((gcs * 100).astype(int), 25, 65)
    exp_by_stratum = {}
    for s in range(25, 66):
        v = obs[ok & (strata == s)]
        if len(v) >= 30:
            exp_by_stratum[s] = np.median(v)
    expected = np.array([exp_by_stratum.get(s, np.nan) for s in strata])
    valid = ~np.isnan(expected) & (obs > 0.3 * med) & (obs < 2.0 * med)
    ratio = obs[valid] / expected[valid]
    vchrom = chroms[valid]
    vbin = bins_[valid]
    print(f"[fit] {valid.sum()} usable bins; GC strata {sorted(exp_by_stratum)}"[:120])

    print("\nchrom   n_bins  ratio   ci95_lo ci95_hi implied_f   verdict")
    rows = []
    for c in CHROMS:
        r = ratio[vchrom == c]
        if len(r) < 20:
            continue
        m = float(np.median(r))
        boots = [np.median(rng.choice(r, len(r))) for _ in range(400)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        exp0 = 0.5 if c in ("X", "Y") else 1.0
        f = 2 * (m - exp0) / exp0
        sig = "" if lo <= exp0 <= hi else ("GAIN?" if m > exp0 else "LOSS?")
        rows.append((c, len(r), m, lo, hi, f, sig))
        print(f"{c:>5} {len(r):8d} {m:7.4f} {lo:7.4f} {hi:7.4f} {f:+9.3f}   {sig}")

    with open(ROOT / "results" / "mosaic_scan2.tsv", "w") as out:
        out.write("chrom\tn_bins\tratio\tci_lo\tci_hi\timplied_f\tverdict\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # segmental check: 10Mb windows vs own chromosome
    print("\n[segments >3 MAD from own chromosome, 10Mb windows]")
    n_seg = 0
    for c in CHROMS:
        sel = vchrom == c
        if sel.sum() < 50:
            continue
        cmed = np.median(ratio[sel])
        cmad = np.median(np.abs(ratio[sel] - cmed)) * 1.4826 or 1e-9
        segs = defaultdict(list)
        for b, r in zip(vbin[sel], ratio[sel]):
            segs[b // 100].append(r)
        for s, v in sorted(segs.items()):
            if len(v) >= 30:
                z = (np.median(v) - cmed) / (cmad / np.sqrt(len(v)))
                if abs(z) > 3 and abs(np.median(v) / cmed - 1) > 0.04:
                    print(f"  chr{c}:{s*10}-{s*10+10}Mb ratio={np.median(v):.3f} "
                          f"(chrom {cmed:.3f}) z={z:.1f} n={len(v)}")
                    n_seg += 1
    if n_seg == 0:
        print("  none")
    print("\nwrote results/mosaic_scan2.tsv")


if __name__ == "__main__":
    main()
