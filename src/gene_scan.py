#!/usr/bin/env python3
"""gene_scan: extract & annotate all proband variants in one or more genes.

Usage:
    python3 gene_scan.py BUB1B                    # one gene
    python3 gene_scan.py BUB1B CEP57 TRIP13       # several genes
    python3 gene_scan.py --panel panels/sac.txt   # gene list file (one symbol per line)
    python3 gene_scan.py BUB1B --all              # include non-PASS variants
    python3 gene_scan.py BUB1B --tsv out.tsv      # also write TSV

Output columns: gene, chr:pos ref>alt, GT, AD, consequence, HGVSc/p, gnomAD AF, ClinVar, flags.
Interesting rows (rare/novel + protein-affecting, or ClinVar P/LP) are marked with '**'.
"""
import argparse
import sys
from pathlib import Path

from common import gene_coords, vcf_variants_in_region, vep_region, RESULTS_DIR

PROTEIN_AFFECTING = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant",
    "splice_region_variant", "protein_altering_variant", "transcript_ablation",
}
RARE_AF = 0.001  # 0.1%


def pick_csq(vep_rec):
    """Prefer the canonical transcript consequence."""
    tcs = vep_rec.get("transcript_consequences") or []
    for tc in tcs:
        if tc.get("canonical"):
            return tc
    return tcs[0] if tcs else {}


def gnomad_af(vep_rec, ref, alt):
    # VEP keys frequencies by *trimmed* allele for indels ('-' for deletions)
    keys = [alt]
    if len(ref) > len(alt) and ref.startswith(alt):
        keys.append("-")
    if len(alt) > len(ref) and alt.startswith(ref):
        keys.append(alt[len(ref):])
    best = None
    for cv in vep_rec.get("colocated_variants", []):
        freqs = cv.get("frequencies", {})
        frs = [freqs[k] for k in keys if k in freqs]
        if not frs and len(freqs) == 1:
            frs = list(freqs.values())
        for fr in frs:
            for k in ("gnomadg", "af"):
                if k in fr:
                    best = fr[k] if best is None else max(best, fr[k])
    return best


def clinvar_sig(vep_rec):
    sigs = set()
    for cv in vep_rec.get("colocated_variants", []):
        sigs.update(cv.get("clin_sig") or [])
    return ",".join(sorted(sigs))


def rsid(vep_rec):
    for cv in vep_rec.get("colocated_variants", []):
        if str(cv.get("id", "")).startswith("rs"):
            return cv["id"]
    return ""


def scan(symbols, pass_only=True):
    rows = []
    for sym in symbols:
        try:
            g = gene_coords(sym)
        except Exception as e:
            print(f"[warn] lookup failed for {sym}: {e}", file=sys.stderr)
            continue
        if not g:
            print(f"[warn] gene not found in Ensembl: {sym}", file=sys.stderr)
            continue
        region = f"{g['chrom']}:{g['start']}-{g['end']}"
        variants = vcf_variants_in_region(region, pass_only=pass_only)
        print(f"[scan] {sym}: {len(variants)} variants", file=sys.stderr)
        if not variants:
            rows.append(dict(gene=sym, note="no variants in region"))
            continue
        vep_in = [f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ." for v in variants]
        try:
            anns = vep_region(vep_in)
        except Exception as e:
            print(f"[warn] VEP failed for {sym}: {e}", file=sys.stderr)
            rows.append(dict(gene=sym, note=f"VEP failed: {e}"))
            continue
        ann_by_key = {}
        for a in anns:
            p = a["input"].split()
            ann_by_key[(p[0], int(p[1]), p[3], p[4])] = a
        for v in variants:
            a = ann_by_key.get((v["chrom"], v["pos"], v["ref"], v["alt"]), {})
            tc = pick_csq(a)
            af = gnomad_af(a, v["ref"], v["alt"])
            csq = a.get("most_severe_consequence", "?")
            clin = clinvar_sig(a)
            clin_set = set(clin.split(",")) if clin else set()
            interesting = (
                (clin_set & {"pathogenic", "likely_pathogenic"}
                 and (af is None or af < 0.05))
                or (csq in PROTEIN_AFFECTING and (af is None or af < RARE_AF))
            )
            rows.append(dict(
                gene=tc.get("gene_symbol") or sym, chrom=v["chrom"], pos=v["pos"],
                ref=v["ref"], alt=v["alt"], gt=v["gt"], ad=v["ad"], filter=v["filter"],
                consequence=csq,
                hgvsc=(tc.get("hgvsc") or "").split(":")[-1],
                hgvsp=(tc.get("hgvsp") or "").split(":")[-1],
                gnomad_af=af, rsid=rsid(a), clinvar=clin, interesting=interesting,
            ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="*", help="gene symbols")
    ap.add_argument("--panel", help="file with one gene symbol per line (# comments ok)")
    ap.add_argument("--all", action="store_true", help="include non-PASS variants")
    ap.add_argument("--tsv", help="write TSV to this path")
    args = ap.parse_args()

    symbols = list(args.genes)
    if args.panel:
        for line in Path(args.panel).read_text().splitlines():
            line = line.split("#")[0].strip()
            if line:
                symbols.append(line)
    if not symbols:
        ap.error("no genes given")

    rows = scan(symbols, pass_only=not args.all)

    hdr = ["", "gene", "variant", "GT", "AD", "consequence", "HGVSc", "HGVSp",
           "gnomAD_AF", "rsid", "ClinVar"]
    lines = ["\t".join(hdr)]
    n_interesting = 0
    for r in rows:
        if "note" in r:
            lines.append(f"\t{r['gene']}\t({r['note']})")
            continue
        mark = "**" if r["interesting"] else ""
        n_interesting += bool(r["interesting"])
        af = "novel" if r["gnomad_af"] is None else f"{r['gnomad_af']:.3g}"
        lines.append("\t".join([mark, r["gene"], f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}",
                                r["gt"], r["ad"], r["consequence"], r["hgvsc"], r["hgvsp"],
                                af, r["rsid"], r["clinvar"] or "-"]))
    print("\n".join(lines))
    print(f"\n{sum(1 for r in rows if 'note' not in r)} variants across {len(symbols)} genes; "
          f"{n_interesting} flagged interesting (**)", file=sys.stderr)

    if args.tsv:
        Path(args.tsv).write_text("\n".join(lines) + "\n")
        print(f"written: {args.tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
