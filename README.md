# KAIBoost

**Boosting with agent nodes: a residual-driven ensemble for rare-disease variant prioritisation.**

First application: [Rare Disease, Real Kid — MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), Track 1.

---

## The idea

Gradient boosting works because each new tree is fitted to what the ensemble so far
gets *wrong* — the residual — and its output is added to a running total. KAIBoost keeps
that skeleton and generalises the node: **a node may be any agent or tool call** — an
annotation lookup, a database query, a model inference, an LLM reasoning step, a whole
sub-pipeline. "K AI" = K AIs as weak learners.

Two design rules keep this a boosting ensemble rather than "a pipeline with LLMs in it":

1. **Every node is fitted to the residual.** A node's input must explicitly include what the
   current ensemble is uncertain or in conflict about — not merely the raw data again. If a
   node would produce the same output whether or not the earlier nodes had run, it is not a
   boosting node.
2. **Every node's output is additive on a shared scale.** Nodes reason freely inside, but must
   emit a number on one common evidence scale (log-likelihood-ratio-style evidence points;
   ACMG scoring is a domain-native instance of exactly this). Free-form reasoning inside,
   fixed arithmetic at the interface.

And one structural rule, borrowed from model trees (M5):

3. **Agents are leaves and feature generators — never split nodes.** Internal splits stay cheap
   deterministic predicates over a feature table, so the structure is reproducible, cacheable and
   fittable. When routing needs judgement, the judgement is first *materialised as a feature
   column*, then split by an ordinary predicate. Agents don't route; agents' outputs route.

## Why it matters here

A rare-disease case is a residual problem by nature. A panel scan explains most of the signal
cheaply; what remains is a handful of ambiguous, conflicting, or uncovered hypotheses — exactly
the thing a second learner should be fitted to. The expensive nodes (genome-wide re-ranking,
splice-effect inference, read-level re-alignment) are spent only where the ensemble is still
uncertain, and each returns evidence points that add into one score.

That score is submitted as **EPCR** (Estimated Probability of Causal Relationship), the
challenge's own unit.

## What this ensemble did

Applied to a single proband's whole-genome VCF, blind to the answer:

| Node | Role | Result |
|---|---|---|
| Panel scan | strong prior learner | Candidate compound-heterozygous pair in *BUB1B* |
| Exomiser (genome-wide) | independent residual learner | Same gene ranked **1 / 399**, same variant pair selected |
| Local SpliceAI (PyTorch port) | resolves the splice hypothesis | Both intronic alternatives scored **0.00** — mechanism excluded |
| Locus structural check | resolves the "large deletion" hypothesis | No low-depth run, no heterozygosity desert, no split reads |
| Mosaic-aneuploidy scan | tests the disease mechanism directly | **Negative**, with a quantified detection bound |
| Adversarial node | attacks the ensemble's own output | Killed a spurious rank-3 compound het (one repeat-region indel split across two records) |

Two independent nodes converging on the same variant pair is what moved the lead hypothesis
from EPCR 0.60 to 0.72; the splice node's clean negatives are what freed the probability mass
that made room for it.

**The method was calibrated before it was trusted.** Ten ClinVar pathogenic variants, selected by
a rule fixed in advance, were spiked into the real VCF one at a time: **10/10 ranked first** with
the matching phenotype, **0/10** under phenotype ablation. The ranking is driven by the phenotype
match, not by variant properties alone.

See [`report/track1_report.md`](report/track1_report.md) for the full evidence chain, the
validation, and the limitations.

## Repository layout

```
src/
  common.py, gene_scan.py, rare_filter.py, variant_report.py   panel-scan nodes
  make_submission.py                                           EPCR assembly -> submission CSV
  phenotype.py                                                 loads HPO terms from local config
  spike_in.py                                                  10-case blind calibration harness
  spliceai_torch.py, splice_node.py                            local SpliceAI (no TensorFlow)
  aneuploidy_scan.py                                           VCF-only mosaic scan
  sample_pass.sh, stream_filter.py                             streaming FASTQ re-alignment
  mosaic_scan2.py, mosaic_scan3.py                             GC / GC x mappability models
  control_pass.sh, control_pass2.sh, case_vs_control.py        matched normal-control benchmark
config/
  phenotype.example.txt                                        template (real phenotype not committed)
report/
  track1_report.md
```

## Reproducing

Requires `bcftools`, `minimap2`, `samtools`, Python 3.12 with `numpy`, `torch`, `pyfaidx`,
`h5py`, `spliceai` (weights only — the models run under the PyTorch port, not TensorFlow),
and Exomiser 14.0.0 with the 2406 hg38 + phenotype bundles.

```bash
cp config/phenotype.example.txt config/phenotype.txt   # then fill in the proband's HPO terms
python3 src/rare_filter.py                             # panel scan
python3 src/spike_in.py                                # calibration: 10 blind spike-ins
python3 src/splice_node.py                             # local SpliceAI over the locus
bash    src/sample_pass.sh                             # streaming re-alignment -> 100 kb bins
python3 src/mosaic_scan3.py                            # GC x mappability, leave-one-chromosome-out
bash    src/control_pass2.sh && python3 src/case_vs_control.py   # normal-control benchmark
```

### One local note that cost an afternoon

TensorFlow aborts at import on this machine (`libc++abi: mutex lock failed`), so the official
SpliceAI `.h5` weights are executed by walking the stored Keras graph in PyTorch
(`src/spliceai_torch.py`). The port is validated before use: on reference sequence it must
recover the annotated splice sites of the gene under study — it scores them at a median
probability of 0.99, and `src/splice_node.py` aborts if that self-test fails.

## Data policy

**This repository contains no patient data, and no file derived from it.**

The challenge dataset is gated under a Data Transfer Agreement: redistribution is prohibited
and all data, including intermediate and derived datasets, is deleted within 30 days of the
hackathon close. The proband's phenotype comes from a gated clinical document and is read at
runtime from `config/phenotype.txt`, which is not committed. `.gitignore` denies every
data-bearing file extension by default.

## Acknowledgement

> This work was made possible through the Hackathon, organized by Sage Bionetworks in
> partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and
> Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are
> deeply grateful to the child and their family who generously contributed their data and their
> story to advance research into this rare disease. We acknowledge their trust in making this
> Hackathon possible.

## Licence

MIT (see [LICENSE](LICENSE)). Hackathon submissions are additionally released under CC-BY 4.0
per the challenge rules.
