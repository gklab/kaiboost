#!/usr/bin/env python3
"""rare_filter: genome-wide funnel — 4.7M PASS variants → rare protein-affecting shortlist.

Pipeline (all intermediates cached under results/):
  1. Build a coding BED (CDS ±8bp splice flank) from the Ensembl GTF (downloaded once, ~50MB).
  2. Restrict proband VCF to coding BED, PASS only            → ~tens of thousands
  3. VEP-annotate those in batches (cached)                   → consequence + gnomAD AF
  4. Keep protein-affecting & rare (AF<0.1%, or <1% if hom)   → shortlist TSV
  5. Score each gene by HPO overlap with the patient's terms  → ranked output

Usage:
    python3 rare_filter.py            # full run (resumes from cached steps)
    python3 rare_filter.py --top 50   # show top 50 by phenotype score
"""
import argparse
import gzip
import subprocess
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

from common import ROOT, CACHE_DIR, RESULTS_DIR, VCF, vep_region, bcftools
from phenotype import load_hpo

GTF_URL = "https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"
G2P_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt"
SPLICE_FLANK = 8
RARE_AF = 0.001
RARE_AF_HOM = 0.01

# Proband HPO terms come from the gated clinical document and must not be
# committed; they are read from config/phenotype.txt (see phenotype.py).
PATIENT_HPO = set(load_hpo())

PROTEIN_AFFECTING = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant",
    "protein_altering_variant", "transcript_ablation",
}


def download(url: str, dest: Path):
    if dest.exists():
        return dest
    print(f"[download] {url}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return dest


def build_coding_bed() -> Path:
    bed = CACHE_DIR / "coding_cds_pm8.bed"
    if bed.exists():
        return bed
    gtf = download(GTF_URL, CACHE_DIR / "ensembl_113.gtf.gz")
    print("[bed] extracting CDS ±8bp from GTF...", file=sys.stderr)
    intervals = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if p[2] != "CDS":
                continue
            # VCF uses no-chr GRCh38 naming, same as Ensembl GTF — no conversion needed
            intervals.append((p[0], max(0, int(p[3]) - 1 - SPLICE_FLANK), int(p[4]) + SPLICE_FLANK))
    # merge
    intervals.sort()
    merged = []
    for c, s, e in intervals:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])
    with open(bed, "w") as f:
        for c, s, e in merged:
            f.write(f"{c}\t{s}\t{e}\n")
    print(f"[bed] {len(merged)} merged intervals", file=sys.stderr)
    return bed


def extract_coding_variants() -> Path:
    out = RESULTS_DIR / "coding_pass.tsv"
    if out.exists():
        return out
    bed = build_coding_bed()
    print("[extract] restricting VCF to coding BED (PASS only)...", file=sys.stderr)
    txt = bcftools(["query", "-R", str(bed), "-i", 'FILTER="PASS"',
                    "-f", r"%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\n", str(VCF)])
    rows, seen = [], set()
    for line in txt.strip().split("\n"):
        c, pos, ref, alts, gt, ad = line.split("\t")
        for alt in alts.split(","):
            key = (c, pos, ref, alt)
            if alt in ("<NON_REF>", "*") or key in seen:
                continue
            seen.add(key)
            rows.append("\t".join([c, pos, ref, alt, gt, ad]))
    out.write_text("\n".join(rows) + "\n")
    print(f"[extract] {len(rows)} coding-region PASS variants", file=sys.stderr)
    return out


def annotate_and_filter() -> Path:
    out = RESULTS_DIR / "rare_shortlist.tsv"
    if out.exists():
        return out
    coding = extract_coding_variants().read_text().strip().split("\n")
    print(f"[vep] annotating {len(coding)} variants in batches of 200 (cached)...", file=sys.stderr)
    shortlist = []
    for i in range(0, len(coding), 200):
        chunk = coding[i:i + 200]
        vs, meta = [], {}
        for line in chunk:
            c, pos, ref, alt, gt, ad = line.split("\t")
            vs.append(f"{c} {pos} . {ref} {alt} . . .")
            meta[(c, int(pos), ref, alt)] = (gt, ad)
        try:
            anns = vep_region(vs)
        except Exception as e:
            print(f"[vep] batch {i // 200} failed: {e} — skipping", file=sys.stderr)
            continue
        for a in anns:
            p = a["input"].split()
            key = (p[0], int(p[1]), p[3], p[4])
            gt, ad = meta.get(key, ("?", "?"))
            csq = a.get("most_severe_consequence", "?")
            if csq not in PROTEIN_AFFECTING:
                continue
            ref_a, alt_a = p[3], p[4]
            keys = [alt_a]
            if len(ref_a) > len(alt_a) and ref_a.startswith(alt_a):
                keys.append("-")
            if len(alt_a) > len(ref_a) and alt_a.startswith(ref_a):
                keys.append(alt_a[len(ref_a):])
            af = None
            for cv in a.get("colocated_variants", []):
                freqs = cv.get("frequencies", {})
                frs = [freqs[k] for k in keys if k in freqs]
                if not frs and len(freqs) == 1:
                    frs = list(freqs.values())
                for fr in frs:
                    for k in ("gnomadg", "af"):
                        if k in fr:
                            af = fr[k] if af is None else max(af, fr[k])
            clin = ",".join(sorted({s for cv in a.get("colocated_variants", [])
                                    for s in (cv.get("clin_sig") or [])}))
            hom = gt in ("1/1", "1|1")
            cutoff = RARE_AF_HOM if hom else RARE_AF
            if af is not None and af >= cutoff:
                continue
            tcs = a.get("transcript_consequences") or [{}]
            tc = next((t for t in tcs if t.get("canonical")), tcs[0])
            shortlist.append([p[0], p[1], p[3], p[4], gt, ad, csq,
                              tc.get("gene_symbol", "?"),
                              (tc.get("hgvsc") or "").split(":")[-1],
                              (tc.get("hgvsp") or "").split(":")[-1],
                              "novel" if af is None else f"{af:.3g}", clin or "-"])
        done = min(i + 200, len(coding))
        if done % 2000 < 200 or done == len(coding):
            print(f"[vep] {done}/{len(coding)} annotated, shortlist={len(shortlist)}",
                  file=sys.stderr)
    hdr = ["chrom", "pos", "ref", "alt", "gt", "ad", "consequence", "gene",
           "hgvsc", "hgvsp", "gnomad_af", "clinvar"]
    out.write_text("\t".join(hdr) + "\n" + "\n".join("\t".join(r) for r in shortlist) + "\n")
    print(f"[filter] shortlist written: {out} ({len(shortlist)} variants)", file=sys.stderr)
    return out


def load_gene_hpo() -> dict[str, set]:
    g2p = download(G2P_URL, CACHE_DIR / "genes_to_phenotype.txt")
    gene_hpo = defaultdict(set)
    with open(g2p) as f:
        header = f.readline().rstrip("\n").split("\t")
        gi = header.index("gene_symbol") if "gene_symbol" in header else 1
        hi = header.index("hpo_id") if "hpo_id" in header else 2
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(gi, hi):
                gene_hpo[p[gi]].add(p[hi])
    return gene_hpo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    shortlist_path = annotate_and_filter()
    rows = [line.split("\t") for line in
            shortlist_path.read_text().strip().split("\n")[1:]]

    gene_hpo = load_gene_hpo()
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r[7]].append(r)

    scored = []
    for gene, vs in by_gene.items():
        overlap = PATIENT_HPO & gene_hpo.get(gene, set())
        n_hits = len(vs)
        has_lof = any(v[6] in ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
                               "splice_donor_variant", "start_lost") for v in vs)
        has_hom = any(v[4] in ("1/1", "1|1") for v in vs)
        biallelic_possible = n_hits >= 2 or has_hom
        score = len(overlap) * 10 + biallelic_possible * 3 + has_lof * 2 + min(n_hits, 3)
        scored.append((score, gene, len(overlap), n_hits, biallelic_possible, has_lof, vs))
    scored.sort(reverse=True)

    print(f"\n=== Top {args.top} genes (phenotype-weighted) — "
          f"{len(rows)} rare protein-affecting variants in {len(by_gene)} genes ===\n")
    print(f"{'score':>5} {'gene':12} {'HPO∩':>4} {'nvar':>4} {'biallelic?':>10} {'LoF?':>4}")
    for score, gene, ov, n, bi, lof, vs in scored[:args.top]:
        print(f"{score:>5} {gene:12} {ov:>4} {n:>4} {str(bi):>10} {str(lof):>4}")
        for v in vs:
            print(f"        {v[0]}:{v[1]} {v[2]}>{v[3]} GT={v[4]} {v[6]} {v[9] or v[8]} "
                  f"AF={v[10]} ClinVar={v[11]}")
    ranked_out = RESULTS_DIR / "ranked_genes.tsv"
    with open(ranked_out, "w") as f:
        f.write("score\tgene\thpo_overlap\tn_variants\tbiallelic\thas_lof\n")
        for score, gene, ov, n, bi, lof, _ in scored:
            f.write(f"{score}\t{gene}\t{ov}\t{n}\t{bi}\t{lof}\n")
    print(f"\nfull ranking: {ranked_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
