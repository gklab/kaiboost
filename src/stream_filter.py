#!/usr/bin/env python3
"""Streaming SAM consumer for the FASTQ sample pass.

Reads SAM on stdin (minimap2 output), never writes a BAM:
  - counts primary MAPQ>=20 reads into 100 kb bins per contig;
  - saves raw SAM lines landing in the BUB1B locus (15:40.05-40.28 Mb);
  - prints a one-line summary to stderr at EOF.

Usage: ... | python3 stream_filter.py <out_prefix>
Writes <out_prefix>.bins.tsv and <out_prefix>.locus.sam (+ header).
"""
import sys
from collections import defaultdict

BIN = 100_000
LOCUS = ("15", 40_050_000, 40_280_000)
MAIN = {str(c) for c in range(1, 23)} | {"X", "Y"}

prefix = sys.argv[1]
bins = defaultdict(int)
locus = open(prefix + ".locus.sam", "w")
n_tot = n_used = n_locus = 0

for line in sys.stdin:
    if line.startswith("@"):
        locus.write(line)
        continue
    try:
        f = line.split("\t", 6)
        flag = int(f[1])
        rname, pos, mapq = f[2], int(f[3]), int(f[4])
    except (IndexError, ValueError):
        continue
    n_tot += 1
    if flag & 0x4:              # unmapped
        continue
    if rname.startswith("chr"):
        rname = rname[3:]
    if rname.startswith("chr"):
        rname = rname[3:]
    # locus file keeps supplementary/secondary too (split reads = SV evidence)
    if rname == LOCUS[0] and LOCUS[1] <= pos <= LOCUS[2]:
        locus.write(line)
        n_locus += 1
    if flag & 0x900 or mapq < 20 or rname not in MAIN:
        continue
    n_used += 1
    bins[(rname, pos // BIN)] += 1

locus.close()
with open(prefix + ".bins.tsv", "w") as out:
    out.write("chrom\tbin\tcount\n")
    for (c, b), n in sorted(bins.items()):
        out.write(f"{c}\t{b}\t{n}\n")
print(f"[stream_filter {prefix}] total={n_tot} used(MQ20 primary)={n_used} "
      f"locus={n_locus}", file=sys.stderr)
