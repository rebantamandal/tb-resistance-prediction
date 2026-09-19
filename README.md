# TB Resistance Prediction

**Predict whether an antibiotic will work against a tuberculosis infection, by
reading the bacterium's DNA.**

Five trained models — rifampicin, isoniazid, ethambutol, pyrazinamide and
streptomycin — plus everything used to build them.

```bash
python predict_resistance.py --variants your_samples.csv --out answers.csv --trust-local-models
```

```text
sample       Rifampicin   score   Isoniazid    score
ERR1034590   Resistant    0.91    Resistant    0.98
ERR047009    Susceptible  0.17    Susceptible  0.16
```

> ### ⚠️ Research use only
>
> These models are **not medically approved**. The scores are **not
> probabilities** — 0.80 does not mean an 80% chance. Never use this output to
> decide anyone's treatment. Growing the bacteria in a laboratory remains the
> real test.

> **New to this, or not from a biology background?**
> **[`EXPLAINED.md`](EXPLAINED.md)** covers the whole project in plain English —
> what it does, how to read the output, and what the results mean.

---

## The problem it addresses

To find out whether an antibiotic works against someone's TB infection, the
standard method is to grow the bacteria in a laboratory alongside the drug and
see whether it survives. **That takes 2 to 8 weeks**, because TB grows slowly.

Reading the DNA takes **1 to 2 days**.

In between, a patient is treated on a guess. This project tests whether the DNA
can give a useful answer in the meantime.

---

## Getting started

Python 3.11 or newer.

```bash
git clone https://github.com/rebantamandal/tb-resistance-prediction.git
cd tb-resistance-prediction

python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Windows
# .venv/bin/python -m pip install -r requirements.txt     # macOS / Linux
```

The trained models are included, so you can predict straight away. **No large
download needed.**

### Try it now

A ready-made example ships with the repo:

```bash
.venv/Scripts/python predict_resistance.py \
  --variants examples/sample_isolates.csv \
  --out test.csv \
  --trust-local-models
```

```text
sample       features  Rifampicin   score   Isoniazid    score
ERR1034590   30        Resistant    0.91    Resistant    0.98
ERR1034591   30        Resistant    0.85    Resistant    0.89
ERR047002    41        Susceptible  0.25    Resistant    0.74
ERR047009     5        Susceptible  0.17    Susceptible  0.16
```

Those four samples have known laboratory results and the models never trained on
them. **All four rifampicin predictions are correct.** See
[`examples/`](examples/).

---

## What you feed it

One CSV. Five columns. **One row per DNA change, per sample.**

```csv
SAMPLE,CHROM,POS,REF,ALT
Patient1,NC_000962.3,761155,C,T
Patient1,NC_000962.3,2155168,C,G
```

| Column | Meaning |
|---|---|
| `SAMPLE` | Your name for the sample |
| `CHROM` | Which reference the positions were measured against |
| `POS` | Where the change is — a position number along the DNA |
| `REF` | What the letter should be |
| `ALT` | What it actually is |

You don't write this by hand; it comes out of the DNA sequencing process.

**Rows are not samples.** One sample takes thousands of rows, because it has
thousands of DNA changes. The example file has 6,781 rows and just 4 samples.

> ⚠️ **Positions must be measured against the reference `NC_000962.3`.** A
> position number only means something relative to a particular reference. With
> the wrong one, the tool recognises nothing and calls everything "Susceptible."
> The `known_features_present` column tells you if this happened — if it's `0`,
> that's the problem.

## What you get back

| Column | Meaning |
|---|---|
| `isolate_id` | Your sample name |
| `known_features_present` | How many recognisable DNA changes were found |
| `<Drug>_call` | **Resistant** (drug will probably fail) or **Susceptible** |
| `<Drug>_score` | Confidence, 0 to 1. At or above 0.50 becomes Resistant. |
| `warning` | Fills in when nothing was recognised |

Scores between **0.4 and 0.6** mean the tool is genuinely unsure — those are the
ones to send for laboratory testing.

---

## How accurate is it?

| Drug | Samples learned from | Familiar hospitals | **New hospital** | New country |
|---|---:|---:|---:|---|
| Rifampicin | 2,625 | 0.958 | **0.908** | 0.69–0.76 |
| Streptomycin | 1,573 | 0.948 | **0.897** | 0.36–0.86 |
| Isoniazid | 2,603 | 0.947 | **0.855** | 0.77–0.94 |
| Ethambutol | 2,408 | 0.913 | **0.765** | 0.63–0.84 |
| Pyrazinamide | 1,839 | 0.911 | **0.782** | 0.63–0.84 |

These are **AUC** scores: show the tool one resistant and one non-resistant
sample — how often does it pick correctly? 0.5 is a coin flip, 1.0 is perfect.

**Plan around the bold column.** Here's why the columns differ:

- **Familiar hospitals** — tested on samples from collections it also trained
  on. This is the figure usually published, and the most flattering: the model
  can partly succeed by recognising the *collection* rather than the disease.
- **New hospital** — tested on collections it has never seen. This costs
  **0.05 to 0.15 across all five drugs**, every time, with no overlap in the
  margins of error. **This is the project's main finding**, and it's what
  realistic use looks like.
- **New country** — a whole country removed from training. Noisy and *not*
  uniformly worse. Isoniazid scores 0.936 on a country it never saw. But
  streptomycin scores **0.359 on South Africa — worse than guessing**, a real
  warning that the tool can fail badly on an unfamiliar population.

Full numbers with margins of error: [`RESULTS.md`](RESULTS.md). Why the columns
differ, explained simply: [`EXPLAINED.md`](EXPLAINED.md).

### It learned real biology

Ask the rifampicin model which DNA position mattered most, and it points to
**the exact mutation medicine has known causes rifampicin resistance for
decades** — about twice as important as anything else. Nobody told it where to
look; it found that from raw position numbers.

### It's been verified

- Re-scoring 150 samples the model was tested on during training: **150 of 150
  identical decisions**, biggest score difference 0.00005 (rounding).
- A fresh download of this repository runs the predictor and all 69 software
  tests with nothing extra to install, producing identical numbers.

---

## How it was built

```text
all_variants.csv          15 million DNA changes, 9,842 samples, NO answers
        |
        |-- prepare_variants.py       find the 814 DNA positions that matter
        |                             (out of 758,768 present)
        |
        |-- fetch_sample_metadata.py  collect sample details from a public database
        |
        |-- fetch_phenotypes.py       recover real laboratory answers from
        |-- build_phenotype_table.py  three public medical databases
        |                             -> 2,664 samples with real answers
        |
        |-- amr.py train              train the models
```

**The source dataset advertised answers it didn't contain.**
[FORUM-TB on Kaggle](https://www.kaggle.com/datasets/nanzhen/forum-tb) describes
"4 drug resistance labels" in a "247 KB" file; what's actually published is
555 MB of DNA changes with no answers at all. The answers here were recovered
independently from BV-BRC, CRyPTIC and NCBI Pathogen Detection — laboratory
results only, never other software's predictions.

---

## What's in this repository

**Start here**

| File | What it's for |
|---|---|
| [`EXPLAINED.md`](EXPLAINED.md) | The whole project in plain English |
| [`USAGE.md`](USAGE.md) | Day-to-day reference: commands, formats, troubleshooting |
| [`RESULTS.md`](RESULTS.md) | Every number, with margins of error |
| [`examples/`](examples/) | A file you can run immediately |

**The tools**

| File | What it does |
|---|---|
| `predict_resistance.py` | **The predictor.** DNA in, answers out. |
| `prepare_variants.py` | Turns raw DNA changes into a table models can use |
| `fetch_sample_metadata.py` | Collects sample details from the ENA public database |
| `fetch_phenotypes.py` | Collects laboratory answers from BV-BRC |
| `build_phenotype_table.py` | Merges answers from three databases, tracking sources |
| `run_full_study.py` | Trains and tests a drug all three ways |
| `geographic_sweep.py` | Tests each country separately |
| `analyze_runs.py` | Calculates margins of error |
| `amr.py` | The underlying model trainer |
| `app.py` | Optional browser interface, if you'd rather not use the command line |

**The data that ships with it**

| Path | What it is |
|---|---|
| `models/*_genomic_baseline/` | The five trained models |
| `data/phenotypes.csv` | Recovered laboratory answers, 15 drugs, with sources |
| `data/sample_metadata.csv` | Sample details for all 9,842 samples |
| `data/variant_panel/feature_dictionary.csv` | The 814 DNA positions the models use |
| `resources/` | The reference tables the pipeline needs |

Working folders are rebuilt by the commands in [`RESULTS.md`](RESULTS.md) §8
rather than stored here, to keep the download small.

---

## Adding another drug

Laboratory answers already exist for several more drugs in `data/phenotypes.csv`
— amikacin (1,163 samples), capreomycin (936), kanamycin (719) and others. No
new data needed:

```bash
python run_full_study.py --drugs AMIKACIN
```

Then add a short entry to `DEFAULT_REGISTRY` in `predict_resistance.py`, copying
an existing one, and it appears in the output alongside the rest.

---

## Known limitations

- **Scores aren't probabilities.** Use them as rankings.
- **The training samples are unusual** — about two thirds carry a major
  resistance mutation, far above a normal population. Accuracy on ordinary
  patients is untested.
- **It can fail badly on an unfamiliar population.** Streptomycin scored below
  chance on South African samples.
- **Positions must use the `NC_000962.3` reference** or nothing is recognised.
- **Only 27% of samples had laboratory answers**, which limits how precise the
  measurements can be.
- **Patient details mostly don't exist.** Age, previous treatment and hospital
  stay aren't recorded publicly for these samples, and weren't invented.
- **Location couldn't be evaluated.** Every collection in this data comes from
  exactly one country, so "which country" and "which hospital" are the same
  information. See [`RESULTS.md`](RESULTS.md) §4.

---

## Data sources and licence

| Source | Used for | Terms |
|---|---|---|
| [FORUM-TB](https://www.kaggle.com/datasets/nanzhen/forum-tb) | DNA changes | CC BY-SA 4.0 |
| [ENA Portal](https://www.ebi.ac.uk/ena/portal/api/) | Sample details | Open |
| [BV-BRC](https://www.bv-brc.org/) | Laboratory answers | Public |
| [CRyPTIC](https://ftp.ebi.ac.uk/pub/databases/cryptic/) | Laboratory answers | See their [data paper](https://doi.org/10.1101/2021.09.14.460274) |
| [NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) | Laboratory answers | Public domain |
| [WHO mutation catalogue](https://www.who.int/publications/i/item/9789240082410) | Reference for which DNA regions matter | — |

Code is MIT licensed. The trained models and DNA-derived data are built from a
CC BY-SA 4.0 source and are shared under those same terms. See
[`LICENSE`](LICENSE).

## Citation

```bibtex
@software{mandal_tb_resistance_prediction,
  author = {Mandal, Rebanta},
  title  = {TB Resistance Prediction: predicting M. tuberculosis drug
            resistance from whole-genome variants},
  url    = {https://github.com/rebantamandal/tb-resistance-prediction},
  year   = {2026}
}
```

Please cite the underlying data sources above as well.
