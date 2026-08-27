#!/usr/bin/env python3
"""Spike-in validation of the variant-ranking method.

For N known ClinVar-pathogenic variants (deterministic selection, no cherry-picking):
  1. spike the variant into the real proband VCF (hom for AR diseases, het for AD),
  2. build a phenopacket from the disease's HPO profile (arm A)
     and one from the MVA proband's unrelated HPO terms (arm B, phenotype ablation),
  3. run Exomiser on both arms,
  4. record the rank of the spiked gene / variant.

Selection rule (fixed BEFORE looking at results): iterate chromosomes 1..22 in order,
take the first variant on each chromosome satisfying ALL of:
  - CLNSIG == Pathogenic (exact), CLNREVSTAT contains 'multiple_submitters'
  - SNV (len(ref)==len(alt)==1), single ALT
  - GENEINFO names exactly one gene; that gene has a disease in genes_to_phenotype
    with an AR/AD inheritance term and >=4 other HPO terms
  - position absent from the proband VCF
  - gene not used already
Stop at N cases. All selected cases are reported, including failures.
"""
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # mva-analysis
EXO = ROOT.parent / "exomiser"
PROBAND_VCF = ROOT.parent / "mva-data" / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
CLINVAR = ROOT / "cache" / "clinvar.vcf.gz"
G2P = ROOT / "cache" / "genes_to_phenotype.txt"
WORK = ROOT / "validation" / "work"
OUT_TSV = ROOT / "validation" / "spike_in_results.tsv"
JAVA = "/opt/homebrew/opt/openjdk/bin/java"
N_CASES = 10

AR, AD = "HP:0000007", "HP:0000006"
INHERITANCE_BRANCH = {AR, AD, "HP:0001417", "HP:0001419", "HP:0001423", "HP:0001450",
                      "HP:0001428", "HP:0032113", "HP:0034345", "HP:0000005"}
sys.path.insert(0, str(ROOT / "tools"))
from phenotype import load_hpo  # noqa: E402

# Arm B (phenotype ablation) uses the proband's own, unrelated phenotype. It is
# read from the uncommitted config/phenotype.txt, never hardcoded, because the
# clinical document it comes from is gated and may not be redistributed.
MVA_HPO = load_hpo()


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))}\n"
                           f"STDERR: {r.stderr[-800:]}\nSTDOUT: {r.stdout[-1500:]}")
    return r.stdout


def load_disease_profiles():
    """gene -> (disease_id, moi, [hpo terms])  choosing the disease with most terms."""
    by_gene_disease = defaultdict(lambda: defaultdict(set))
    with open(G2P) as f:
        header = f.readline().rstrip("\n").split("\t")
        gi = header.index("gene_symbol")
        hi = header.index("hpo_id")
        di = header.index("disease_id")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(gi, hi, di):
                by_gene_disease[p[gi]][p[di]].add(p[hi])
    profiles = {}
    for gene, diseases in by_gene_disease.items():
        best = None
        for did, terms in diseases.items():
            moi = "AR" if AR in terms else ("AD" if AD in terms else None)
            pheno = sorted(t for t in terms if t not in INHERITANCE_BRANCH)
            if moi and len(pheno) >= 4:
                cand = (len(pheno), did, moi, pheno[:10])
                if best is None or cand[0] > best[0]:
                    best = cand
        if best:
            profiles[gene] = (best[1], best[2], best[3])
    return profiles


def in_proband(chrom, pos):
    out = run(["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}", "-f", "%POS\n",
               str(PROBAND_VCF)])
    return str(pos) in out.split()


def select_cases(profiles):
    cases, used_genes = [], set()
    for chrom in [str(c) for c in range(1, 23)]:
        if len(cases) >= N_CASES:
            break
        stream = subprocess.Popen(
            ["bcftools", "query", "-r", chrom,
             "-i", 'INFO/CLNSIG="Pathogenic" && INFO/CLNREVSTAT~"multiple_submitters"',
             "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENEINFO\n", str(CLINVAR)],
            stdout=subprocess.PIPE, text=True)
        try:
            for line in stream.stdout:
                c, pos, ref, alt, geneinfo = line.rstrip("\n").split("\t")
                if len(ref) != 1 or len(alt) != 1 or alt not in "ACGT":
                    continue
                genes = [g.split(":")[0] for g in geneinfo.split("|")]
                if len(genes) != 1 or genes[0] in used_genes or genes[0] not in profiles:
                    continue
                gene = genes[0]
                if in_proband(c, pos):
                    continue
                did, moi, pheno = profiles[gene]
                cases.append(dict(chrom=c, pos=int(pos), ref=ref, alt=alt, gene=gene,
                                  disease=did, moi=moi, hpo=pheno))
                used_genes.add(gene)
                break
        finally:
            stream.terminate()
    return cases


def make_spiked_vcf(case, dest_dir):
    gt = "1/1:0,45:45:99" if case["moi"] == "AR" else "0/1:22,23:45:99"
    row = (f"{case['chrom']}\t{case['pos']}\t.\t{case['ref']}\t{case['alt']}\t1000\tPASS\t"
           f"DP=45;MQ=60;QD=25\tGT:AD:DP:GQ\t{gt}\n")
    spike = dest_dir / "spike.vcf"
    with open(spike, "w") as f:
        f.write(run(["bcftools", "view", "-h", str(PROBAND_VCF)]))
        f.write(row)
    run(["bcftools", "view", "-Oz", "-o", f"{spike}.gz", str(spike)])
    run(["bcftools", "index", "-t", f"{spike}.gz"])
    merged = dest_dir / "merged.vcf.gz"
    run(["bcftools", "concat", "-a", "-Oz", "-o", str(merged), str(PROBAND_VCF),
         f"{spike}.gz"])
    run(["bcftools", "index", "-t", str(merged)])
    return merged


def make_phenopacket(case_id, hpo_terms, path):
    # subject id MUST match the VCF sample name or Exomiser rejects the job
    lines = ["---", f"id: {case_id}", "subject:", "  id: WGS_EX2312012", "  sex: MALE",
             "phenotypicFeatures:"]
    for t in hpo_terms:
        lines += [f"  - type:", f"      id: {t}", f"      label: {t}"]
    lines += ["metaData:", "  created: '2026-08-27T00:00:00.000Z'", "  createdBy: spikein",
              "  resources:", "    - id: hp", "      name: human phenotype ontology",
              "      url: http://purl.obolibrary.org/obo/hp.owl",
              "      version: hp/releases/2024-06-01", "      namespacePrefix: HP",
              "      iriPrefix: 'http://purl.obolibrary.org/obo/HP_'",
              "  phenopacketSchemaVersion: 1.0"]
    path.write_text("\n".join(lines) + "\n")


def run_exomiser(vcf, sample_yml, outdir, name):
    run([JAVA, "-Xmx16g", "-jar", "exomiser-cli-14.0.0.jar",
         "--sample", str(sample_yml), "--vcf", str(vcf), "--assembly", "hg38",
         "--preset", "exome", "--output-directory", str(outdir),
         "--output-filename", name, "--output-format", "TSV_GENE,TSV_VARIANT"],
        cwd=EXO / "exomiser-cli-14.0.0")


def gene_rank(genes_tsv, gene):
    best = None
    with open(genes_tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            if row["GENE_SYMBOL"] == gene:
                r = int(row["#RANK"])
                best = r if best is None else min(best, r)
    return best


def variant_contributing(variants_tsv, case):
    with open(variants_tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            if (row["CONTIG"].lstrip("chr") == case["chrom"]
                    and row["START"] == str(case["pos"]) and row["REF"] == case["ref"]
                    and row["ALT"] == case["alt"] and row["CONTRIBUTING_VARIANT"] == "1"):
                return True
    return False


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    print("[1/3] loading disease profiles...", flush=True)
    profiles = load_disease_profiles()
    print(f"      {len(profiles)} genes with usable disease profiles", flush=True)
    print("[2/3] selecting cases (deterministic rule)...", flush=True)
    cases = select_cases(profiles)
    for c in cases:
        print(f"      {c['gene']:10} {c['chrom']}:{c['pos']} {c['ref']}>{c['alt']} "
              f"{c['moi']} {c['disease']} ({len(c['hpo'])} HPO)", flush=True)

    results = []
    print("[3/3] running spike-in cases...", flush=True)
    for i, case in enumerate(cases):
        cdir = WORK / f"case{i:02d}_{case['gene']}"
        cdir.mkdir(parents=True, exist_ok=True)
        try:
            merged = make_spiked_vcf(case, cdir)
            row = dict(case=f"case{i:02d}", **{k: case[k] for k in
                       ("gene", "chrom", "pos", "ref", "alt", "moi", "disease")})
            for arm, hpo in (("A_true_pheno", case["hpo"]), ("B_mva_pheno", MVA_HPO)):
                yml = cdir / f"{arm}.yml"
                make_phenopacket(f"case{i:02d}", hpo, yml)
                run_exomiser(merged, yml, cdir, arm)
                g = gene_rank(cdir / f"{arm}.genes.tsv", case["gene"])
                contrib = variant_contributing(cdir / f"{arm}.variants.tsv", case)
                row[f"{arm}_gene_rank"] = g
                row[f"{arm}_variant_contributing"] = contrib
                print(f"      {row['case']} {case['gene']:10} {arm}: gene_rank={g} "
                      f"contributing={contrib}", flush=True)
            results.append(row)
        except Exception as e:
            print(f"      {case['gene']} FAILED: {e}", flush=True)
            results.append(dict(case=f"case{i:02d}", gene=case["gene"], error=str(e)[:200]))
        finally:
            for f in cdir.glob("*.vcf.gz*"):
                f.unlink()
            (cdir / "spike.vcf").unlink(missing_ok=True)

    keys = sorted({k for r in results for k in r}, key=lambda k: (k != "case", k))
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    ok = [r for r in results if "error" not in r]
    ranks = [r["A_true_pheno_gene_rank"] for r in ok if r.get("A_true_pheno_gene_rank")]
    if ranks:
        t1 = sum(1 for r in ranks if r == 1)
        t5 = sum(1 for r in ranks if r <= 5)
        t10 = sum(1 for r in ranks if r <= 10)
        print(f"\n=== ARM A (true phenotype): top1 {t1}/{len(ranks)}  "
              f"top5 {t5}/{len(ranks)}  top10 {t10}/{len(ranks)}", flush=True)
    ranks_b = [r["B_mva_pheno_gene_rank"] for r in ok if r.get("B_mva_pheno_gene_rank")]
    if ranks_b:
        t1b = sum(1 for r in ranks_b if r == 1)
        print(f"=== ARM B (wrong phenotype): top1 {t1b}/{len(ranks_b)}  "
              f"ranks: {ranks_b}", flush=True)
    print(f"results: {OUT_TSV}", flush=True)


if __name__ == "__main__":
    main()
