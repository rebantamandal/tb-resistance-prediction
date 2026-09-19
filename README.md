# TB Resistance Prediction

Predict *Mycobacterium tuberculosis* drug resistance from whole-genome variant
calls. Four trained Random Forest models — rifampicin, isoniazid, ethambutol and
pyrazinamide — plus the full pipeline that built them, from a 15-million-row
variant table to a one-command predictor.

```bash
python predict_resistance.py --variants your_isolates.csv --out predictions.csv --trust-local-models
```

```text
isolate_id   Rifampicin_call  Rifampicin_score  Isoniazid_call  Isoniazid_score
ERR038736    Resistant        0.5360            Resistant       0.8330
SRR998857    Resistant        0.8552            Resistant       0.9418
ERR038266    Susceptible      0.0592            Susceptible     0.1320
```

> ### ⚠️ Research use only
>
> These models are **not clinically validated** and their scores are
> **uncalibrated**. A score of 0.80 is not an 80% probability that an infection
> is resistant. Do not use any output here to select, start, stop or change
> antibiotic treatment. Culture-based drug susceptibility testing remains the
> reference standard.

> **New to this, or not from a biology background?** Read
> **[`EXPLAINED.md`](EXPLAINED.md)** — the whole project in plain English, with
> no jargon: what the tool does, how to read its output, and what the results
> mean.

---

## Why this exists

Tuberculosis drug susceptibility testing takes **2–8 weeks** because TB grows
slowly. Sequencing takes **1–2 days**. In between, patients are treated on a
guess — and if the guess is wrong they stay sick and infectious.

Resistance in TB is driven by mutations in a small, known set of genes (`rpoB`
for rifampicin, `katG` for isoniazid, and so on), which is what makes prediction
from a genome tractable at all.

---

## Quick start

Python 3.11 or newer.

```bash
git clone https://github.com/rebantamandal/tb-resistance-prediction.git
cd tb-resistance-prediction

python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Windows
# .venv/bin/python -m pip install -r requirements.txt     # macOS / Linux
```

The four trained models are included in the clone, so you can predict
immediately — no training and no 555 MB download required:

```bash
.venv/Scripts/python predict_resistance.py \
  --variants your_isolates.csv \
  --out predictions.csv \
  --trust-local-models
```

`--trust-local-models` is required every time and deliberately so: model files
execute code when loaded, so you confirm each run that these are files you trust.

### Try it right now

The repository ships four real isolates with known laboratory results:

```bash
.venv/Scripts/python predict_resistance.py   --variants examples/sample_isolates.csv   --out my_predictions.csv   --trust-local-models
```

```text
isolate_id   known_features  Rifampicin   score   Isoniazid    score
ERR1034590   30              Resistant    0.9062  Resistant    0.9797
ERR1034591   30              Resistant    0.8491  Resistant    0.8948
ERR047002    41              Susceptible  0.2473  Resistant    0.7360
ERR047009     5              Susceptible  0.1746  Susceptible  0.1584
```

All four were in the rifampicin model's held-out test set, and all four
rifampicin calls match the real laboratory result. See
[`examples/`](examples/) for details.

---

## Input format

One CSV, five columns, one row per variant call per isolate. Any number of
isolates per file.

```csv
SAMPLE,CHROM,POS,REF,ALT
ERR9001,NC_000962.3,761155,C,T
ERR9001,NC_000962.3,2155168,C,G
ERR9002,NC_000962.3,1673425,C,T
```

| Column | Meaning |
|---|---|
| `SAMPLE` | Your isolate identifier. Anything unique; echoed back unchanged. |
| `CHROM` | Reference contig. Not used, but see the warning below. |
| `POS` | 1-based position on the reference. |
| `REF` | Reference allele. |
| `ALT` | Observed allele. |

> **Positions must be called against H37Rv (`NC_000962.3`).** A different
> reference produces an all-zero feature profile that scores like a susceptible
> isolate. The tool flags any isolate carrying no recognised feature — treat that
> flag as "check your coordinates", never as a susceptible result.

## Output

One row per isolate, with a score and a call for every drug:

| Column | Meaning |
|---|---|
| `isolate_id` | Your `SAMPLE` value |
| `variant_calls_supplied` | Rows you provided for this isolate |
| `known_features_present` | How many of the model's 814 features it carries |
| `<Drug>_score` | Uncalibrated resistance score, 0–1 |
| `<Drug>_call` | `Resistant` at or above the threshold (default 0.50) |
| `warning` | Set when no known resistance feature was detected |

A `.report.json` is written alongside, recording the threshold, feature count,
flagged isolates, and the measured performance of every model used.

---

## How well does it work?

Held-out ROC-AUC with 95% bootstrap intervals. **Same models, same features —
only the train/test split changes.**

| Drug | Labels | Grouped by isolate | Grouped by study | Held-out country |
|---|---:|---:|---:|---|
| Rifampicin | 2,625 | 0.958 | **0.908** | 0.69–0.76 |
| Isoniazid | 2,603 | 0.947 | **0.855** | 0.77–0.94 |
| Ethambutol | 2,408 | 0.913 | **0.765** | 0.63–0.84 |
| Pyrazinamide | 1,839 | 0.911 | **0.782** | 0.63–0.84 |

**Use the bold column as your working expectation.** It is what the model scores
on isolates from a study it never trained on — the closest match to routine use.

### What the columns mean

- **Grouped by isolate** — the protocol usually reported in the literature.
  Isolates from one study land on both sides of the split, so the model can
  partly recognise the cohort rather than the biology.
- **Grouped by study** — no study straddles the split. This costs **0.05–0.15
  AUC across all four drugs**, with non-overlapping intervals every time. This
  is the project's most robust finding.
- **Held-out country** — an entire country excluded from training. Results here
  are **mixed and noisy**, not a uniform collapse. Rifampicin drops on all three
  countries tested (Canada 0.686 [0.565, 0.799]). Isoniazid does the opposite,
  scoring 0.936 on held-out Canada against 0.855 within-cohort. Six of fourteen
  country results beat their study-grouped figure, and several intervals span
  0.4 or more.

Full per-country results with intervals are in [`RESULTS.md`](RESULTS.md) §2.

### It learns real biology

Permutation importance on the rifampicin model ranks **rpoB Ser450Leu** — the
dominant rifampicin resistance mutation worldwide — far above every other
feature, at roughly twice the next. The pipeline recovered it from raw genomic
coordinates without being told what to look for.

### Verified against the training pipeline

Scoring 150 isolates from the rifampicin model's own held-out set through the
standalone predictor:

- maximum score difference from the recorded values: **5 × 10⁻⁵** (output rounding)
- resistant/susceptible calls: **150 of 150 identical**
- against real laboratory labels: accuracy **0.847**, ROC-AUC **0.947**

A fresh `git clone` of this repository was verified to run the predictor and the
full test suite with no extra downloads, producing byte-identical scores.

---

## How it was built

```text
all_variants.csv                    15,057,917 variant calls, 9,842 isolates, no labels
        │
        ├─ prepare_variants.py      long → wide; 43 resistance loci + 200 bp promoter
        │                           margins; drop variants in <20 isolates
        │                           → 9,842 × 814 binary feature matrix
        │
        ├─ fetch_sample_metadata.py real country / city / specimen from the ENA API
        │
        ├─ fetch_phenotypes.py      laboratory DST from BV-BRC (lab evidence only)
        ├─ build_phenotype_table.py + CRyPTIC + NCBI Pathogen Detection
        │                           → 2,664 isolates with real labels
        │
        └─ amr.py train             Random Forest, 300 trees, group-disjoint split
                                    → the four models in models/
```

**The source dataset ships no labels.** [FORUM-TB on Kaggle](https://www.kaggle.com/datasets/nanzhen/forum-tb)
advertises "4 drug resistance labels" and "247 KB compressed", but the published
file is 555 MB of unlabelled variant calls. Labels were recovered independently
from three public sources of laboratory susceptibility testing, joined by
sequencing run accession:

| Source | Isolates |
|---|---:|
| BV-BRC `genome_amr` (laboratory evidence only) | 2,408 |
| CRyPTIC reuse table | 257 |
| NCBI Pathogen Detection `AST_phenotypes` | 79 |

Computational predictions were excluded at the query — training on another
model's output would be circular. Where BV-BRC and NCBI cover the same
isolate-drug pair, they **agree on 323 of 323 results**, so label noise is not
what limits the numbers above.

---

## Repository layout

| Path | What it is |
|---|---|
| `predict_resistance.py` | **The predictor.** Variants in, per-drug predictions out. |
| `prepare_variants.py` | Builds the isolate × variant feature matrix from raw calls. |
| `fetch_sample_metadata.py` | Real country, city and specimen from the ENA portal API. |
| `fetch_phenotypes.py` | Laboratory DST from BV-BRC. |
| `build_phenotype_table.py` | Merges the three label sources, with provenance and agreement. |
| `run_full_study.py` | Runs the whole evaluation ladder for any drug. |
| `geographic_sweep.py` | Holds out each country in turn. |
| `analyze_runs.py` | Paired bootstrap confidence intervals on saved runs. |
| `amr.py` | The underlying Random Forest trainer and CLI. |
| `app.py` | Optional local browser interface for training and prediction. |
| `models/*_genomic_baseline/` | The four released models, with their reports. |
| `resources/tb_resistance_loci.csv` | 43 H37Rv resistance loci with drug associations. |
| `data/phenotypes.csv` | Recovered laboratory labels, 15 drugs, with per-label source. |
| `data/sample_metadata.csv` | ENA metadata for all 9,842 accessions. |
| `data/variant_panel/feature_dictionary.csv` | Every feature's position, alleles and gene. |

Documentation: [`EXPLAINED.md`](EXPLAINED.md) (plain-English guide, start here) ·
[`USAGE.md`](USAGE.md) (day-to-day use) ·
[`RESULTS.md`](RESULTS.md) (what was measured and how) ·
[`DATA_PIPELINE.md`](DATA_PIPELINE.md) (how the data was assembled) ·
[`TESTING.md`](TESTING.md).

---

## Reproducing from scratch

Download `all_variants.csv` from [FORUM-TB](https://www.kaggle.com/datasets/nanzhen/forum-tb)
into the project root, then:

```bash
# 1. Feature matrix (~5 min, two streaming passes over 555 MB)
python prepare_variants.py --variants all_variants.csv \
  --out-dir data/variant_panel --min-prevalence 20

# 2. Metadata and labels (public APIs, no account needed)
python fetch_sample_metadata.py \
  --accessions data/variant_panel/isolate_variant_matrix.csv \
  --out data/sample_metadata.csv
python fetch_phenotypes.py --out data/phenotype_sources/bvbrc_phenotypes.csv \
  --accessions data/variant_panel/isolate_variant_matrix.csv
python build_phenotype_table.py \
  --bvbrc data/phenotype_sources/bvbrc_phenotypes.csv \
  --accessions data/variant_panel/isolate_variant_matrix.csv \
  --out data/phenotypes.csv

# 3. Train and evaluate any drug
python run_full_study.py --drugs RIFAMPICIN ISONIAZID ETHAMBUTOL PYRAZINAMIDE

# 4. Confidence intervals
python analyze_runs.py --runs models/rif_genomic_baseline
```

Labels already exist for eleven more drugs in `data/phenotypes.csv`
(streptomycin 1,573; amikacin 1,163; capreomycin 936; and others). Add one with
`run_full_study.py --drugs STREPTOMYCIN`, then register it in
`DEFAULT_REGISTRY` in `predict_resistance.py`.

### Tests

```bash
python -m unittest discover -s tests
```

69 tests, covering feature-name reproducibility across processes, locus
annotation margins, label handling, metadata normalisation, phenotype merging
and the predictor's matrix construction.

---

## Limitations

- **Not calibrated.** Scores rank isolates; they are not probabilities.
- **Trained on an unusual population.** 65% of isolates carry the katG isoniazid
  mutation. Resistance prevalence here is far above community TB, so precision
  will not transfer to a screening setting unchanged.
- **Drops sharply across geography.** ~0.91 within a familiar population,
  ~0.69–0.84 on an unseen country. Validate locally before trusting it anywhere new.
- **H37Rv coordinates only.** Other references silently yield empty profiles.
- **No call-quality filtering.** The source carries no depth or genotype-quality
  fields, so every variant call is taken at face value.
- **No patient identifiers** exist in the source, so repeat isolates from one
  patient may cross a split, flattering every figure slightly.
- **27% label coverage** — 2,664 of 9,842 isolates. More labels would tighten
  every interval.
- **Clinical features are mostly unavailable.** Age, prior antibiotic exposure
  and hospitalisation duration are not recorded for these accessions in any
  public source, and were not invented. Only infection site was obtainable.
- **Geography and clinical features could not be evaluated.** Every study in
  this cohort sits in exactly one country (Cramér's V between study and country
  = 1.000), so "which country" and "which cohort" are the same variable. Group
  the split by isolate and country becomes a cohort lookup (+0.017 AUC); group
  it by study and country is an unseen category carrying nothing (−0.039 AUC).
  Both intervals exclude zero and they disagree in sign, because each measures a
  different artefact. The released models are genomic-only for that reason — not
  because geography is irrelevant to resistance, which remains untested here.
  See [`RESULTS.md`](RESULTS.md) §3.

---

## Data sources and attribution

| Source | Used for | Terms |
|---|---|---|
| [FORUM-TB](https://www.kaggle.com/datasets/nanzhen/forum-tb) | Variant calls | CC BY-SA 4.0 |
| [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/) | Sample metadata | Open |
| [BV-BRC](https://www.bv-brc.org/) | Laboratory DST | Public |
| [CRyPTIC](https://ftp.ebi.ac.uk/pub/databases/cryptic/) | Laboratory DST | See their [data compendium](https://doi.org/10.1101/2021.09.14.460274) |
| [NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) | Laboratory DST | Public domain |
| [WHO mutation catalogue, 2nd ed.](https://www.who.int/publications/i/item/9789240082410) | Locus annotation reference | — |

The source variant data is CC BY-SA 4.0. The feature matrix, feature dictionary
and trained models are derived from it and are shared under the **same
CC BY-SA 4.0 terms**. The code in this repository is MIT licensed — see
[`LICENSE`](LICENSE).

## Citation

```bibtex
@software{mandal_tb_resistance_prediction,
  author = {Mandal, Rebanta},
  title  = {TB Resistance Prediction: Random Forest models for
            M. tuberculosis drug resistance from whole-genome variants},
  url    = {https://github.com/rebantamandal/tb-resistance-prediction},
  year   = {2026}
}
```

Please also cite the underlying data sources above, particularly the CRyPTIC
data compendium and the FORUM-TB dataset.
