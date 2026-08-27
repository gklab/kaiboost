#!/usr/bin/env python3
"""KAIBoost aneuploidy node: detect mosaic whole-chromosome gains/losses from
the proband VCF alone (no BAM), using two aggregate signals per chromosome:

  1. depth ratio: median DP at variant sites vs autosome-wide median
     (mosaic trisomy fraction f -> ratio (2+f)/2; monosomy -> (2-f)/2);
  2. het allele-balance overdispersion: trisomy splits het AB into branches
     (1+f)/(2+f) and 1/(2+f); excess sd beyond binomial implies
     f = 4d/(1-2d) where d is the branch half-offset.

Male proband: chrX (haploid) must show depth ratio ~0.5 — internal control.
Also emits a 10 Mb-bin depth table to flag large segmental CNVs, and a
focused report on the BUB1B locus (15:40.10-40.23 Mb) for a het deletion
in trans with the truncating allele.
"""
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
VCF = ROOT.parent / "mva-data" / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
OUT = ROOT / "results" / "aneuploidy_scan.tsv"
BINS_OUT = ROOT / "results" / "depth_bins_10mb.tsv"
CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]
DP_LO, DP_HI = 25, 80
BIN = 10_000_000
BUB1B = ("15", 40_100_000, 40_230_000)


def stream(args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    yield from p.stdout
    p.wait()


def main():
    het_ab = defaultdict(list)      # chrom -> AB of PASS het SNVs (DP window)
    het_dp = defaultdict(list)
    all_dp_bins = defaultdict(list)  # (chrom, bin) -> DP of all PASS sites
    locus = []                       # (pos, gt, ab, dp) inside BUB1B locus

    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD\t%DP]\n",
           "-i", 'FILTER="PASS"', str(VCF)]
    n = 0
    for line in stream(cmd):
        f = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, gt, ad, dp = f[0], int(f[1]), f[2], f[3], f[4], f[5], f[6]
        if chrom not in ("15",) and chrom not in CHROMS:
            continue
        n += 1
        try:
            dpv = int(dp)
        except ValueError:
            continue
        all_dp_bins[(chrom, pos // BIN)].append(dpv)
        is_snv = len(ref) == 1 and len(alt) == 1 and alt in "ACGT"
        het = gt in ("0/1", "0|1", "1/0", "1|0")
        ab = None
        if het and is_snv and DP_LO <= dpv <= DP_HI:
            parts = ad.split(",")
            if len(parts) == 2:
                r, a = int(parts[0]), int(parts[1])
                if r + a >= DP_LO:
                    ab = a / (r + a)
                    het_ab[chrom].append(ab)
                    het_dp[chrom].append(dpv)
        if chrom == BUB1B[0] and BUB1B[1] <= pos <= BUB1B[2]:
            locus.append((pos, gt, -1.0 if ab is None else ab, dpv))
        if n % 1_000_000 == 0:
            print(f"  {n/1e6:.0f}M sites...", flush=True)

    # per-chromosome table
    auto_dp = np.concatenate([het_dp[c] for c in CHROMS[:22] if het_dp[c]])
    genome_med_dp = np.median(auto_dp)
    rows = []
    for c in CHROMS:
        if not het_ab[c]:
            rows.append([c, 0] + ["."] * 7)
            continue
        ab = np.array(het_ab[c])
        dp = np.array(het_dp[c], dtype=float)
        med_dp = float(np.median(dp))
        ratio = med_dp / genome_med_dp
        sd = float(ab.std())
        binom_sd = float(np.sqrt(np.mean(0.25 / dp)))
        excess = float(np.sqrt(max(sd**2 - binom_sd**2, 0.0)))
        f_tri_ab = 4 * excess / (1 - 2 * excess) if excess < 0.5 else float("nan")
        f_dp = 2 * (ratio - 1)
        rows.append([c, len(ab), round(float(ab.mean()), 4), round(sd, 4),
                     round(binom_sd, 4), round(excess, 4), round(f_tri_ab, 3),
                     round(med_dp, 1), round(ratio, 3), round(f_dp, 3)])

    hdr = ["chrom", "n_het", "mean_ab", "sd_ab", "binom_sd", "excess_sd",
           "f_implied_ab", "med_dp", "dp_ratio", "f_implied_dp"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # 10Mb depth bins, flag |dev| > 20% vs own chromosome median
    with open(BINS_OUT, "w") as fh:
        fh.write("chrom\tbin_start_mb\tn\tmed_dp\tvs_chrom\n")
        flagged = []
        for c in CHROMS:
            cbins = {b: v for (cc, b), v in all_dp_bins.items() if cc == c}
            if not cbins:
                continue
            cmed = np.median([np.median(v) for v in cbins.values()])
            for b in sorted(cbins):
                m = float(np.median(cbins[b]))
                dev = m / cmed
                fh.write(f"{c}\t{b*10}\t{len(cbins[b])}\t{m:.1f}\t{dev:.3f}\n")
                if abs(dev - 1) > 0.20 and len(cbins[b]) > 200:
                    flagged.append((c, b * 10, m, dev))

    print("\n=== per-chromosome (sorted by |dp_ratio-1|):")
    print("\t".join(hdr))
    for r in sorted(rows, key=lambda r: -abs((r[8] if r[8] != "." else 1) - 1)):
        print("\t".join(str(x) for x in r))
    print(f"\n=== 10Mb bins deviating >20% from own chromosome: {len(flagged)}")
    for c, mb, m, dev in flagged[:20]:
        print(f"  chr{c}:{mb}Mb med_dp={m:.0f} ratio={dev:.2f}")

    # BUB1B locus: het deletion in trans would show halved DP + het desert
    ldp = np.array([d for _, _, _, d in locus], dtype=float)
    lhet = [x for x in locus if x[2] >= 0]
    print(f"\n=== BUB1B locus 15:40.10-40.23Mb: {len(locus)} PASS sites, "
          f"{len(lhet)} usable hets")
    print(f"  med_dp={np.median(ldp):.1f} (genome {genome_med_dp:.1f}, "
          f"ratio {np.median(ldp)/genome_med_dp:.2f})")
    if lhet:
        abs_ = np.array([x[2] for x in lhet])
        print(f"  het AB mean={abs_.mean():.3f} sd={abs_.std():.3f} n={len(abs_)}")
    print(f"\nwrote {OUT} and {BINS_OUT}")


if __name__ == "__main__":
    main()
