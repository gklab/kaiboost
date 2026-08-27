#!/usr/bin/env python3
"""variant_report: deep-dive dossier for a single variant.

Usage:
    python3 variant_report.py 15:40209701:T:G
    python3 variant_report.py 15:40209701:T:G 15:40220612:T:G   # several at once

Pulls: proband VCF record, VEP annotation (HGVS, consequence, exon), and MyVariant.info
(gnomAD genome+exome AF, ClinVar, CADD, REVEL, AlphaMissense, SpliceAI).
"""
import sys

from common import vcf_variants_in_region, vep_region, myvariant


def dig(d, *path):
    for k in path:
        if isinstance(d, dict):
            d = d.get(k)
        else:
            return None
    return d


def as_scalar(x):
    """MyVariant returns scalar or list depending on transcript count."""
    if isinstance(x, list):
        x = [v for v in x if v is not None]
        return max(x) if x and all(isinstance(v, (int, float)) for v in x) else (x[0] if x else None)
    return x


def report(spec: str):
    chrom, pos, ref, alt = spec.replace("-", ":").split(":")
    chrom = chrom.removeprefix("chr")
    pos = int(pos)
    print(f"\n{'=' * 70}\n VARIANT {chrom}:{pos} {ref}>{alt} (GRCh38)\n{'=' * 70}")

    # 1. proband VCF record
    recs = [v for v in vcf_variants_in_region(f"{chrom}:{pos}-{pos}")
            if v["pos"] == pos and v["ref"] == ref and v["alt"] == alt]
    if recs:
        v = recs[0]
        print(f"[proband]   GT={v['gt']}  AD={v['ad']} (ref,alt reads)  DP={v['dp']}  "
              f"GQ={v['gq']}  FILTER={v['filter']}")
    else:
        print("[proband]   NOT PRESENT in proband VCF")

    # 2. VEP
    a = (vep_region([f"{chrom} {pos} . {ref} {alt} . . ."]) or [{}])[0]
    tc = next((t for t in a.get("transcript_consequences", []) if t.get("canonical")),
              (a.get("transcript_consequences") or [{}])[0])
    print(f"[VEP]       gene={tc.get('gene_symbol', '?')}  "
          f"consequence={a.get('most_severe_consequence', '?')}  "
          f"exon={tc.get('exon') or tc.get('intron') or '-'}")
    print(f"            HGVSc={tc.get('hgvsc', '-')}")
    print(f"            HGVSp={tc.get('hgvsp', '-')}")

    # 3. MyVariant.info
    hgvs = f"chr{chrom}:g.{pos}{ref}>{alt}" if len(ref) == 1 and len(alt) == 1 else None
    if hgvs is None:  # indel HGVS
        if len(ref) > len(alt):  # deletion
            hgvs = f"chr{chrom}:g.{pos + 1}_{pos + len(ref) - 1}del"
        else:  # insertion
            hgvs = f"chr{chrom}:g.{pos}_{pos + 1}ins{alt[1:]}"
    m = (myvariant([hgvs]) or [{}])[0]
    if m.get("notfound"):
        print(f"[MyVariant] not found ({hgvs}) — likely absent from all databases (novel)")
        return
    g_af = as_scalar(dig(m, "gnomad_genome", "af", "af"))
    e_af = as_scalar(dig(m, "gnomad_exome", "af", "af"))
    clin = dig(m, "clinvar", "clinical_significance")
    rcv = dig(m, "clinvar", "rcv")
    if clin is None and rcv:
        rcvs = rcv if isinstance(rcv, list) else [rcv]
        clin = ";".join(sorted({str(r.get("clinical_significance")) for r in rcvs}))
    print(f"[gnomAD]    genome_AF={g_af if g_af is not None else 'absent'}  "
          f"exome_AF={e_af if e_af is not None else 'absent'}")
    print(f"[ClinVar]   {clin or 'not in ClinVar'}")
    print(f"[scores]    CADD={as_scalar(dig(m, 'cadd', 'phred')) or '-'}  "
          f"REVEL={as_scalar(dig(m, 'dbnsfp', 'revel', 'score')) or '-'}  "
          f"AlphaMissense={as_scalar(dig(m, 'dbnsfp', 'alphamissense', 'score')) or '-'}  "
          f"SpliceAI={as_scalar(dig(m, 'dbnsfp', 'spliceai', 'ds_max')) or '-'}")
    print(f"[rsid]      {as_scalar(dig(m, 'dbsnp', 'rsid')) or '-'}")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for spec in sys.argv[1:]:
        report(spec)


if __name__ == "__main__":
    main()
