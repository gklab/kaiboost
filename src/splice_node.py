#!/usr/bin/env python3
"""KAIBoost splice node: SpliceAI delta scores for every PASS variant in the
BUB1B locus, computed fully locally (torch port of the official weights).

Steps:
  0. self-test: reference-sequence acceptor/donor probabilities must peak at
     the annotated BUB1B exon boundaries (validates the weight conversion);
  1. score mva-analysis/splice/bub1b_locus.vcf at -D 4999, raw and masked;
  2. write mva-analysis/splice/bub1b_spliceai.tsv sorted by max raw DS.
"""
import sys
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from spliceai_torch import TorchSpliceAI

# stub keras before importing spliceai.utils (TF is unusable in this env)
keras_mod = types.ModuleType("keras")
keras_models = types.ModuleType("keras.models")
keras_models.load_model = TorchSpliceAI
keras_mod.models = keras_models
sys.modules["keras"] = keras_mod
sys.modules["keras.models"] = keras_models

import spliceai.utils as U  # noqa: E402


def one_hot_encode(seq):  # np.fromstring was removed in numpy 2.x
    m = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0],
                    [0, 0, 1, 0], [0, 0, 0, 1]])
    seq = seq.upper().replace("A", "\x01").replace("C", "\x02") \
             .replace("G", "\x03").replace("T", "\x04").replace("N", "\x00")
    return m[np.frombuffer(seq.encode("latin1"), np.int8) % 5]


U.one_hot_encode = one_hot_encode

FASTA = ROOT / "cache" / "chr15.fa"
VCF = ROOT / "splice" / "bub1b_locus.vcf"
OUT = ROOT / "splice" / "bub1b_spliceai.tsv"
DIST = 4999


class Record:
    def __init__(self, chrom, pos, ref, alts):
        self.chrom, self.pos, self.ref, self.alts = chrom, pos, ref, alts


def self_test(ann):
    """Reference BUB1B sequence: acceptor/donor must peak at exon boundaries."""
    i = next(k for k, g in enumerate(ann.genes)
             if g == "BUB1B" and ann.chroms[k].lstrip("chr") == "15")
    starts, ends = ann.exon_starts[i], ann.exon_ends[i]
    chrom = U.normalise_chrom("15", list(ann.ref_fasta.keys())[0])
    lo, hi = int(starts[0]) - 5001, int(ends[-1]) + 5000
    seq = ann.ref_fasta[chrom][lo:hi].seq
    y = np.mean([m.predict(one_hot_encode(seq)[None, :]) for m in ann.models],
                axis=0)[0]                      # (L-10000, 3), pos0 = lo+5001
    off = lo + 5001                             # y[i] <-> 1-based position off+i

    def peak(pos1, ch):                         # max prob within +/-2 of site
        i = int(pos1) - off
        return float(y[max(i - 2, 0):i + 3, ch].max())

    acc = [peak(s, 1) for s in starts[1:]]      # internal acceptors
    don = [peak(e, 2) for e in ends[:-1]]       # internal donors
    print(f"[self-test] BUB1B internal splice sites (n={len(acc)}+{len(don)}):")
    print(f"  acceptor probs: median={np.median(acc):.3f} min={min(acc):.3f}")
    print(f"  donor probs:    median={np.median(don):.3f} min={min(don):.3f}")
    if np.median(acc) < 0.5 or np.median(don) < 0.5:
        sys.exit("SELF-TEST FAILED — weight conversion is wrong, aborting")
    print("[self-test] PASSED\n")


def parse_vcf(path):
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        yield Record(f[0], int(f[1]), f[3], tuple(f[4].split(","))), f[9].split(":")[0]


def main():
    print("[load] 5 models (torch)...", flush=True)
    ann = U.Annotator(str(FASTA), "grch38")
    self_test(ann)

    rows = []
    recs = list(parse_vcf(VCF))
    for n, (rec, gt) in enumerate(recs, 1):
        raw = U.get_delta_scores(rec, ann, DIST, 0)
        msk = U.get_delta_scores(rec, ann, DIST, 1)
        for r, m in zip(raw, msk):
            p = r.split("|")
            if p[2] == ".":
                continue
            ds = [float(x) for x in p[2:6]]
            rows.append([rec.chrom, rec.pos, rec.ref, p[0], p[1], gt,
                         max(ds), *p[2:10], *m.split("|")[2:6]])
        if n % 10 == 0:
            print(f"  {n}/{len(recs)} variants scored", flush=True)

    rows.sort(key=lambda r: -r[6])
    hdr = ["chrom", "pos", "ref", "alt", "gene", "gt", "ds_max",
           "ds_ag", "ds_al", "ds_dg", "ds_dl", "dp_ag", "dp_al", "dp_dg", "dp_dl",
           "m_ag", "m_al", "m_dg", "m_dl"]
    with open(OUT, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print(f"\n[done] {len(rows)} scored rows -> {OUT}")
    for r in rows[:8]:
        print(f"  {r[0]}:{r[1]} {r[2]}>{r[3]} {r[4]} gt={r[5]} "
              f"ds_max={r[6]:.2f} (AG {r[7]} AL {r[8]} DG {r[9]} DL {r[10]})")


if __name__ == "__main__":
    main()
