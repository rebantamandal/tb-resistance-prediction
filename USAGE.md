# How to use the predictor

Everything you need to run the tool day to day. No biology background assumed —
where a term is unavoidable, it's explained where it appears.

New to the project entirely? Read [`EXPLAINED.md`](EXPLAINED.md) first.

> **Research use only.** These models are not medically approved and their
> scores are not probabilities. Never use the output to decide someone's
> treatment.

---

## The command

```bat
.venv\Scripts\python predict_resistance.py ^
  --variants your_file.csv ^
  --out answers.csv ^
  --trust-local-models
```

On macOS or Linux use `.venv/bin/python`, and `\` instead of `^` for line breaks.

**Why `--trust-local-models` is required every time:** the saved model files can
run code when they're opened. Typing the flag is you confirming these are your
own files and not something downloaded from a stranger. It's a safety
speed-bump, not a formality.

---

## What goes in

One CSV file. Five columns. **One row per DNA change, per sample.**

```csv
SAMPLE,CHROM,POS,REF,ALT
Patient1,NC_000962.3,761155,C,T
Patient1,NC_000962.3,2155168,C,G
Patient2,NC_000962.3,1673425,C,T
```

| Column | What it is |
|---|---|
| `SAMPLE` | Your name for the sample. Anything unique. Comes back unchanged. |
| `CHROM` | Which reference the positions were measured against. See the warning below. |
| `POS` | *Where* the change is — position number along the DNA. |
| `REF` | What the letter should be. |
| `ALT` | What it actually is. |

You won't write this file yourself. It comes out of the DNA sequencing process.

### Rows are not samples

The most common point of confusion. **A single sample takes up thousands of
rows**, because it has thousands of DNA changes. The example file shipped with
this repo has 6,781 rows and only 4 samples:

```text
ERR047002     2,496 rows
ERR047009       728 rows
ERR1034590    1,774 rows
ERR1034591    1,783 rows
              -----
              6,781 rows  ->  4 samples
```

A file for 100 patients would be roughly 170,000 rows.

### The one thing that silently breaks it

**Positions must be measured against the reference called `NC_000962.3`.**

Position 761,155 only means something relative to a particular reference
document. Measured against a different one, the same number points somewhere
else entirely — and the tool will find nothing it recognises and call everything
"Susceptible."

That's why the output has a `known_features_present` column. If it reads `0`,
you have this problem. A `warning` column fills in to tell you.

---

## What comes out

One row per sample:

| Column | What it means |
|---|---|
| `isolate_id` | Your sample name, echoed back |
| `variant_calls_supplied` | How many rows you gave for this sample |
| `known_features_present` | How many of those the model actually recognised |
| `<Drug>_call` | **Resistant** (drug will probably fail) or **Susceptible** (will probably work) |
| `<Drug>_score` | Confidence, 0 to 1. At or above 0.50 becomes Resistant |
| `warning` | Fills in when nothing was recognised |

A second file, `answers.report.json`, records the settings used and each model's
measured accuracy.

### How to read a score

`0.80` does **not** mean "80% chance this is resistant." The models were never
tuned to produce percentages. Use the number as a ranking — higher means more
confident — and nothing finer.

- **Below 0.4** — fairly confident the drug works
- **0.4 to 0.6** — the tool is genuinely unsure. These are the ones to send for
  laboratory testing.
- **Above 0.6** — fairly confident the drug fails

`--threshold 0.3` catches more resistance at the cost of more false alarms.
`--threshold 0.7` does the reverse. Neither makes the model better; you are
choosing which kind of mistake you would rather make.

---

## The models

| Drug | Samples it learned from | Familiar data | **New hospital** | New country |
|---|---:|---:|---:|---|
| Rifampicin | 2,625 | 0.958 | **0.908** | 0.69-0.76 |
| Isoniazid | 2,603 | 0.947 | **0.855** | 0.77-0.94 |
| Ethambutol | 2,408 | 0.913 | **0.765** | 0.63-0.84 |
| Pyrazinamide | 1,839 | 0.911 | **0.782** | 0.63-0.84 |

These numbers are **AUC**: show the tool one resistant and one non-resistant
sample — how often does it pick correctly? 0.5 is a coin flip, 1.0 is perfect.

**Plan around the bold column.** The three columns differ because of how the
model was tested, not because anything about the model changed:

- **Familiar data** — tested on samples from collections it also trained on.
  This is the number usually published, and the most flattering, because the
  model can partly succeed by recognising the collection rather than the disease.
- **New hospital** — tested on samples from collections it never saw. Costs
  0.05 to 0.15 across all four drugs, every time. **This is realistic use.**
- **New country** — a whole country removed from training. These results are
  noisy and *not* uniformly worse: rifampicin drops on every country tested, but
  isoniazid scores 0.936 on a country it never saw, better than the 0.855 it
  manages on familiar data. Treat any single country number as weak evidence.

[`RESULTS.md`](RESULTS.md) has every figure with its margin of error.
[`EXPLAINED.md`](EXPLAINED.md) explains *why* the columns differ, in plain terms.

All models share one set of 814 recognised positions, so they run together in a
single pass. If you retrain one, retrain them all — the tool refuses to mix
models built from different runs, because their inputs would not line up.

---

## Proof it works

Running the tool on 150 samples the rifampicin model had been tested on during
training, then comparing:

- biggest difference in score: **0.00005** — just decimal rounding
- Resistant/Susceptible decisions: **150 out of 150 identical**
- against the real laboratory results: **84.7% correct**

A fresh download of this repository was also checked: it runs the predictor and
the full test suite with nothing extra to install, producing identical numbers.

Run the software tests yourself:

```bat
.venv\Scripts\python -m unittest discover -s tests
```

69 tests. They check the software behaves correctly — not that the predictions
are medically accurate.

---

## Adding another drug

Laboratory results already exist in `data/phenotypes.csv` for several drugs
beyond the ones shipped, including amikacin (1,163 samples), capreomycin (936)
and kanamycin (719). No new data needed.

```bat
.venv\Scripts\python run_full_study.py --drugs AMIKACIN
```

That trains the model and measures it three ways. Then add a few lines to
`DEFAULT_REGISTRY` near the top of `predict_resistance.py`, copying the pattern
of an existing entry, and it appears in the output alongside the others.

---

## When something goes wrong

| What you see | What it means |
|---|---|
| Every sample says `Susceptible` and `known_features_present` is 0 | Your positions are not measured against `NC_000962.3`. By far the most common problem. |
| `expects different features from` | Your models came from different training runs. Retrain them together. |
| `Environment differs from the one these models were fitted in` | Your installed library versions differ from the ones used for training. Reinstall from `requirements.txt`, or retrain. |
| `already exists` | The tool never overwrites results. Pick a new output filename. |
| `Pass --trust-local-models` | You left the safety flag off. Add it. |
