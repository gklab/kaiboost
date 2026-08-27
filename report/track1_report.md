# Track 1 — Variant Identification

**Rare Disease, Real Kid: MVA Hackathon 2026**
Method: **KAIBoost** — boosting with agent nodes
Repository: <https://github.com/gklab/kaiboost>

---

## 1. Answer

A compound-heterozygous pair in **BUB1B** (OMIM 257300, mosaic variegated aneuploidy syndrome 1),
submitted at EPCR 0.72:

| | Variant | Consequence | Classification |
|---|---|---|---|
| Allele 1 | `chr15:40209701 T>G` | `c.2210T>G` **p.(Leu737Ter)** | ClinVar P/LP; ACMG **Pathogenic** (PVS1, PM2, PP4, PP5_Strong) |
| Allele 2 | `chr15:40220612 T>G` | `c.3006T>G` **p.(Asn1002Lys)** | Not in ClinVar; ACMG **VUS** (PM2, PP4); MVP 0.85, REVEL 0.47; kinase domain |

*BUB1B* encodes BubR1, the core kinase of the spindle assembly checkpoint. Biallelic loss produces
chromosome mis-segregation, constitutional aneuploidy, growth restriction and childhood cancer
predisposition. Every axis of the proband's supplied HPO profile — the oncological, the growth, the
renal and the reproductive-history terms alike — falls inside that single mechanism, which is why
this gene outranks candidates that explain only one axis.

*(The proband's clinical features are not reproduced here. They come from the challenge's gated
clinical document, which may not be redistributed; judges and organisers hold it already.)*

Three lower-ranked alternative second alleles and four secondary findings were also submitted;
their reasoning is in the submission's `notes` column.

---

## 2. Method: what makes this an ensemble rather than a pipeline

KAIBoost keeps gradient boosting's skeleton and generalises the node: **any agent or tool call may
be a node**. Two rules make it boosting rather than sequencing:

1. **Each node is fitted to the residual.** A node's input must include what the current ensemble
   is uncertain or in conflict about. A node that would return the same thing whether or not the
   earlier nodes ran is not a boosting node.
2. **Each node's output is additive on one scale.** Nodes reason freely inside but must emit a
   number on a shared evidence scale. ACMG scoring is a domain-native instance of exactly this;
   here the running total is the challenge's own EPCR.

One structural rule, borrowed from model trees (M5): **agents are leaves and feature generators,
never split nodes.** Internal splits stay cheap deterministic predicates over a feature table, so
the structure is reproducible and cacheable. Where routing needs judgement, the judgement is first
materialised as a feature column and then split by an ordinary predicate. Agents don't route;
agents' outputs route.

The practical consequence: the expensive nodes (genome-wide re-ranking, splice inference, read-level
re-alignment, a matched control benchmark) were spent **only where the ensemble was still uncertain**,
not uniformly over the genome.

---

## 3. Evidence chain

### Node 1 — panel scan (strong prior learner)

Rare-variant filtering over chromosome-instability and paediatric-cancer gene panels, with
phenotype weighting. Output: the *BUB1B* pair above. Residual left behind: *which* variant is the
second allele, and whether the pair is real or a lookalike.

### Node 2 — Exomiser, genome-wide (independent residual learner)

Exomiser 14.0.0 (hg38, 2406 bundles, exome preset) was run on the full VCF with the proband's HPO
profile — deliberately **not** restricted to any panel, so it could disagree.

- 4,962,048 variants → 1,267 passing frequency/pathogenicity filters → 837 in 399 genes
- **`BUB1B_AR` ranks 1 of 399** (combined 0.9751, phenotype 0.8134, p < 0.0001)
- The two contributing variants it selected are **exactly the pair from Node 1**

Two nodes built on different evidence (panel priors vs. genome-wide phenotype-driven ranking)
converging on the same pair is the single largest reduction of residual in this analysis. EPCR
0.60 → 0.72.

### Node 3 — local SpliceAI (resolves the splice hypothesis)

Two alternative second alleles were intronic, and both rested on an unproven splice mechanism.
Rather than leave them as "pending", the splice question was answered locally.

TensorFlow aborts at import on this machine, so the official SpliceAI `.h5` weights are executed
by walking the stored Keras graph in PyTorch (`src/spliceai_torch.py`). **The port is validated
before use**: on reference sequence it must recover the annotated splice sites of *BUB1B* — 22
internal acceptors and 22 internal donors scored at median probability **0.991 / 0.990**, and the
node aborts if that self-test fails.

All 92 PASS variants in `chr15:40,111,000–40,230,000` were then scored at `-D 4999`:

| Variant | Hypothesis | Δ score (AG/AL/DG/DL) | Outcome |
|---|---|---|---|
| `40192892 C>T` | cryptic splice effect | 0.00 / 0.00 / 0.00 / 0.00 | **excluded** |
| `40216470 A>G` | deep-intronic pseudoexon | 0.00 / 0.00 / 0.00 / 0.00 | **excluded** |
| `40209701 T>G` | (allele 1) | max 0.03 | no splice confound — pure truncation |
| `40220612 T>G` | (allele 2) | max 0.02 | no splice confound — pure missense |

The highest score anywhere in the locus was 0.03, so there is also **no unnoticed cryptic-splice
variant** competing to be the second allele. Two hypotheses' probability mass was released to the
leading one; EPCR for those rows dropped to 0.01.

### Node 4 — locus structural check (resolves the "large deletion" hypothesis)

A heterozygous deletion in *trans* would also explain the phenotype. Across 106 PASS sites in the
locus there is no low-depth run and no heterozygosity desert; the two candidate alleles are
separated by continuous heterozygous calls. At read level, 2,006 reads over the locus contain
**zero supplementary alignments** (split reads are the signature of a breakpoint) and 17/2,006
discordant pairs (0.85%, background). A ≥10 kb heterozygous deletion is excluded.

*Honest weight:* at ~1× locus coverage from the sampled re-alignment, absence of split reads is
weak on its own; it is reported as consistent with — not independent proof alongside — the
depth-based exclusion.

### Node 5 — adversarial node (attacks the ensemble's own output)

Exomiser's rank-3 gene, *FANCD2*, presented an apparently compelling splice compound heterozygote.
It is an artefact: `chr3:10046723 AG>A` and `chr3:10046725 TAAG>T` fall two bases apart inside an
AAG repeat — one alignment event represented as two records, therefore in *cis*, not a *trans*
pair. Both are ClinVar Benign and the Fanconi phenotype is absent. Dismissed, and recorded as a
secondary row so the reasoning is visible rather than silently dropped.

---

## 4. Calibration: the method was tested before it was trusted

A rank-1 result means nothing unless the ranker has been shown to produce rank-1 for the right
reason. Ten ClinVar pathogenic variants were selected by a rule **fixed before any result was
seen** (iterate chromosomes 1–22; first variant per chromosome that is a single-nucleotide
`Pathogenic` call with multiple submitters, in a single gene with an AR/AD disease profile of ≥4
HPO terms, and absent from the proband VCF), spiked one at a time into the real VCF, and run
through the same pipeline in two arms.

| Arm | Phenotype supplied | Genes ranked 1st |
|---|---|---|
| **A** | the spiked disease's own HPO profile | **10 / 10** |
| **B** | the proband's (unrelated) phenotype | **0 / 10** — ranks 28, 13, 3, 18, 3, 16, 23, 19, 7, 19 |

Genes: *AGRN, TPO, TRNT1, PIGG, SDHA, FOXC1, DNAAF5, CLN8, DOCK8, ZMYND11*. Every selected case is
reported, including failures; none were discarded.

Arm B is the important one. It shows the ranking is driven by the **phenotype match**, not by
variant properties alone — the same variant, with the wrong phenotype, falls out of the top 10.

---

## 5. A negative result, reported in full

Mosaic variegated aneuploidy is named for a cell-level phenomenon. Testing for it directly — rather
than inferring it from the genotype — seemed worth the effort, so the disease mechanism itself was
made a node.

**Pass 1 (VCF only).** Per-chromosome depth ratios and heterozygous allele-balance overdispersion.
Flagged chr16/17/19/22. All were GC artefacts.

**Pass 2 (read level).** 13.3M read pairs (a 2.5% uniform sample of all four lanes) were streamed
from the raw FASTQs, aligned with minimap2 to hg38 no-alt, and binned at 100 kb — without ever
storing a BAM. A GC-corrected model flagged chr16/17/18/19/21 at +1.4% to +3.2% implied mosaic
fraction. **But chrX — one copy in this male, so its truth is known — came out 1.6% off.** The
"signal" was the size of the method's own demonstrated error, and every flagged 10 Mb segment sat
in pericentromeric or segmentally-duplicated sequence. Both point at mappability, not copy number.

**Pass 3 (mappability-aware).** Coverage was modelled on GC **and** repeat content (the soft-masked
fraction of the reference, a free RepeatMasker proxy), fitted robustly in log space, and — to stop a
chromosome in its own corner of feature space from explaining away its own deviation — fitted
**leave-one-chromosome-out**. On the truth anchor this works: **chrX lands at 0.4997 against a true
0.5000, an error of 0.02%.**

**The decisive test.** GIAB HG002 — a karyotypically normal male — was put through the identical
pipeline, and the proband was binomially thinned to the control's depth (8 seeds) so dispersion
could be compared fairly. Variegated aneuploidy means *excess chromosome-to-chromosome dispersion*,
so that is the statistic:

| | Proband (depth-matched) | Control (HG002) |
|---|---|---|
| Autosomal chromosome-ratio dispersion | **0.57%** | **0.85%** |
| chrX error against known truth | 0.02% | 2.43% |

The proband's chromosome-level spread is **lower** than a normal sample's (0.68×), and the two
samples' per-chromosome deviations do not correlate (r = −0.37) — i.e. these deviations are
sample-specific technical noise, not a fixed artefact that subtracting a control would remove.

**Conclusion: no mosaic aneuploidy is detectable, at a bound of ~1% coverage shift ≈ 2% of cells
trisomic for a given chromosome.** This is "not detected", not "not present": MVA mosaicism is
classically assayed by karyotype in *cultured lymphocytes*, and the fraction in bulk blood DNA
(largely granulocytes) may sit below this bound. No clinical karyotype was provided with the
dataset. The control is also a different platform (HiSeq 2500 2×148 PCR-free vs. the proband's
NovaSeq), which bounds the noise floor without controlling for platform.

The result did not change the answer. It is reported because a mechanism that was actively tested
and not found is a different epistemic state from one that was never looked for.

---

## 6. Limitations

- **Phase is unconfirmed, and cannot be confirmed from this dataset.** The two alleles are 10.9 kb
  apart. A 300k-pair pilot alignment puts the library's insert size at p50 337 bp, p99 746 bp, with
  **zero pairs ≥2 kb**; the intervening heterozygous gaps are 6.8 kb and 4.1 kb. Read-backed phasing
  is physically impossible here, and no parental samples are available. *In trans* is inferred from
  the recessive disease model and the phenotype fit, not demonstrated.
- **The second allele is a VUS.** p.(Asn1002Lys) is novel, in the kinase domain, with supportive but
  not decisive in-silico support (MVP 0.85, REVEL 0.47). Functional work — a BubR1 kinase or
  checkpoint assay — is the natural next step.
- **One hypothesis remains untested rather than excluded.** The alternative second allele held at
  EPCR 0.20 is `chr15:40114056 A>G` (rs191579534, ~47 kb upstream), which matches the profile of the
  hypomorphic regulatory allele described by Ochiai et al. (2014, PNAS) in *BUB1B* MVA cases with
  only one identifiable coding mutation. There is no local assay for a regulatory allele, so it is
  carried as an open alternative rather than dismissed.
- **Mosaicism below ~2% of cells is not excluded** (Section 5).

---

## 7. Reproducibility and data handling

All code is in <https://github.com/gklab/kaiboost>, with per-node entry points listed in the README.
The analysis used only local resources — a local VCF via `bcftools`, locally downloaded public
reference databases (ClinVar, HPO, Ensembl, GRCh38), local Exomiser bundles, and SpliceAI weights
run locally — with no per-variant queries to external annotation services.

The repository contains **no patient data and no file derived from it**. The proband's phenotype is
read at runtime from an uncommitted local config; `.gitignore` denies every data-bearing file
extension by default. All challenge data and derived intermediates will be deleted within 30 days
of the hackathon close and confirmed by email to the organisers, per the Data Transfer Agreement.

---

## 8. Acknowledgement

> This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership
> with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment
> Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to
> the child and their family who generously contributed their data and their story to advance
> research into this rare disease. We acknowledge their trust in making this Hackathon possible.
