#!/usr/bin/env python3
"""make_submission: build & validate a Track 1 predictions CSV.

The official schema (one row per proposed causal variant or compound-het pair):
    proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2,
    epcr, finding_type, notes
Coordinates GRCh38, chromosomes as 'chr15'. Max 10 rows. epcr in (0,1].

Usage:
    # validate an existing csv
    python3 make_submission.py --validate my_submission.csv

    # build from a simple candidates file (see candidates.example.tsv) and validate
    python3 make_submission.py --build candidates.tsv --out team_approach.csv

Candidates file: TSV with columns
    variant1  variant2(or -)  epcr  finding_type  notes
    e.g.  15:40209701:T:G  15:40220612:T:G  0.9  primary  BUB1B compound het
Every variant is cross-checked against the proband VCF (position+alleles must exist).
"""
import argparse
import csv
import sys
from pathlib import Path

from common import vcf_variants_in_region, PROBAND_ID

COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
           "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
MAIN_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M", "MT"]}


def in_proband_vcf(chrom, pos, ref, alt):
    c = chrom.removeprefix("chr")
    try:
        recs = vcf_variants_in_region(f"{c}:{pos}-{pos}")
    except RuntimeError:
        return False
    return any(v["pos"] == int(pos) and v["ref"] == ref and v["alt"] == alt for v in recs)


def validate(path: Path) -> bool:
    errors, warnings = [], []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != COLUMNS:
            errors.append(f"header mismatch\n  expected: {COLUMNS}\n  got:      {reader.fieldnames}")
            rows = []
        else:
            rows = list(reader)

    if len(rows) > 10:
        errors.append(f"{len(rows)} rows — max 10 candidate rows accepted")
    if not rows:
        errors.append("no data rows")

    prev_epcr = None
    for i, r in enumerate(rows, start=2):  # line numbers incl. header
        tag = f"line {i}"
        if r["proband_id"] != PROBAND_ID:
            warnings.append(f"{tag}: proband_id '{r['proband_id']}' != expected '{PROBAND_ID}'")
        # variant 1 (required)
        for chrom, pos, ref, alt, sfx in [
            (r["chrom_1"], r["pos_1"], r["ref_1"], r["alt_1"], "1"),
            (r["chrom_2"], r["pos_2"], r["ref_2"], r["alt_2"], "2"),
        ]:
            blank = not any([chrom, pos, ref, alt])
            if blank:
                if sfx == "1":
                    errors.append(f"{tag}: variant 1 is required")
                continue
            if not all([chrom, pos, ref, alt]):
                errors.append(f"{tag}: variant {sfx} partially filled")
                continue
            if chrom not in MAIN_CHROMS:
                errors.append(f"{tag}: chrom_{sfx} '{chrom}' must look like 'chr15'")
                continue
            if not pos.isdigit():
                errors.append(f"{tag}: pos_{sfx} '{pos}' not an integer")
                continue
            if not (ref.isalpha() and alt.isalpha()):
                errors.append(f"{tag}: ref/alt_{sfx} must be bases (got '{ref}'>'{alt}')")
                continue
            if not in_proband_vcf(chrom, pos, ref, alt):
                errors.append(f"{tag}: {chrom}:{pos} {ref}>{alt} NOT found in proband VCF "
                              f"(wrong build/coords/alleles?)")
        # epcr
        try:
            e = float(r["epcr"])
            if not (0 < e <= 1):
                errors.append(f"{tag}: epcr {e} outside (0, 1]")
            if prev_epcr is not None and e > prev_epcr:
                warnings.append(f"{tag}: epcr {e} higher than previous row — rows should be "
                                f"ranked best-first")
            prev_epcr = e
        except ValueError:
            errors.append(f"{tag}: epcr '{r['epcr']}' not a float")
        if r["finding_type"] not in ("primary", "secondary"):
            errors.append(f"{tag}: finding_type '{r['finding_type']}' must be primary|secondary")

    for w in warnings:
        print(f"  [warn]  {w}")
    for e in errors:
        print(f"  [ERROR] {e}")
    ok = not errors
    print(f"\n{path}: {'VALID' if ok else 'INVALID'} "
          f"({len(rows)} rows, {len(errors)} errors, {len(warnings)} warnings)")
    return ok


def build(cand_path: Path, out_path: Path):
    rows = []
    for line in cand_path.read_text().splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        v1, v2, epcr, ftype, *notes = line.split("\t")
        c1, p1, r1, a1 = v1.split(":")
        c2 = p2 = r2 = a2 = ""
        if v2 not in ("-", ""):
            c2, p2, r2, a2 = v2.split(":")
        add_chr = lambda c: c if c.startswith("chr") or not c else f"chr{c}"
        rows.append([PROBAND_ID, add_chr(c1), p1, r1, a1, add_chr(c2), p2, r2, a2,
                     epcr, ftype, notes[0] if notes else ""])
    rows.sort(key=lambda r: -float(r[9]))
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        w.writerows(rows)
    print(f"written: {out_path} ({len(rows)} rows)\n")
    validate(out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", metavar="CSV")
    ap.add_argument("--build", metavar="CANDIDATES_TSV")
    ap.add_argument("--out", metavar="CSV", help="output for --build")
    args = ap.parse_args()
    if args.validate:
        sys.exit(0 if validate(Path(args.validate)) else 1)
    elif args.build:
        if not args.out:
            ap.error("--build requires --out")
        build(Path(args.build), Path(args.out))
    else:
        ap.error("use --validate or --build")


if __name__ == "__main__":
    main()
