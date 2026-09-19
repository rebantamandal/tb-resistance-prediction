# Using the resistance predictor

**Research use only.** The models are not clinically validated and the scores are
uncalibrated. Do not use output to start, stop or change antibiotic treatment.

> Looking for a plain-English explanation with no jargon? See
> [`EXPLAINED.md`](EXPLAINED.md).

## The one command

```bat
.venv\Scripts\python predict_resistance.py ^
  --variants your_isolates.csv --out predictions.csv --trust-local-models
```

`--trust-local-models` is required and deliberate: model files execute code when
loaded, so you confirm each time that these were produced locally by you.

## Input

One CSV, five columns, one row per variant call per isolate. Any number of
isolates in one file.

```text
SAMPLE,CHROM,POS,REF,ALT
ERR9001,NC_000962.3,761155,C,T
ERR9001,NC_000962.3,2155168,C,G
ERR9002,NC_000962.3,1673425,C,T
```

| Column | Meaning |
|---|---|
| `SAMPLE` | Your isolate identifier. Anything unique; it is echoed back unchanged. |
| `CHROM` | Reference contig. Ignored, but positions **must** be called against H37Rv `NC_000962.3`. |
| `POS` | 1-based position on that reference. |
| `REF` | Reference allele. |
| `ALT` | Observed allele. |

This is exactly the format of `all_variants.csv`, so anything that produced the
training data produces valid input.

**Positions must come from H37Rv.** A different reference silently yields an
all-zero profile that scores like a susceptible isolate. The tool flags isolates
with no recognised feature for exactly this reason — treat that flag as "check
your coordinates", not as a susceptible result.

## Output

One row per isolate:

| Column | Meaning |
|---|---|
| `isolate_id` | Your `SAMPLE` value. |
| `variant_calls_supplied` | How many rows you gave for this isolate. |
| `known_features_present` | How many of the model's 814 features it carries. |
| `<Drug>_score` | Uncalibrated resistance score, 0 to 1. |
| `<Drug>_call` | `Resistant` if the score is at or above the threshold (default 0.50). |
| `warning` | Set when no known feature was detected. |

A `.report.json` beside it records the threshold, the feature count, how many
isolates were flagged, and the measured performance of each model used.

### Reading a score

`0.80` is **not** an 80% probability that this infection is resistant. The models
are not calibrated. What the score supports is ranking and a threshold decision,
nothing finer. Raising `--threshold` trades sensitivity for specificity; it does
not make the model better.

## The four models

| Drug | Training labels | Grouped by isolate | Grouped by study | Held-out country |
|---|---:|---:|---:|---|
| Rifampicin | 2,625 | 0.958 | **0.908** | 0.69&ndash;0.76 |
| Isoniazid | 2,603 | 0.947 | **0.855** | 0.77&ndash;0.94 |
| Ethambutol | 2,408 | 0.913 | **0.765** | 0.63&ndash;0.84 |
| Pyrazinamide | 1,839 | 0.911 | **0.782** | 0.63&ndash;0.84 |

Held-out ROC-AUC. **Use the bold column as your working expectation** &mdash; it is
what the model scores on isolates from a study it did not train on, which is the
closest match to routine use. The first column is the figure usually published
and is the most optimistic; moving to the second costs 0.05&ndash;0.15 AUC
consistently across all four drugs.

The held-out-country column is a **range across several countries, not a
ceiling**. It is noisy and not uniformly worse: rifampicin drops on every
country tested, while isoniazid scores 0.936 on held-out Canada against 0.855
within-cohort. Several intervals span 0.4 or more. See `RESULTS.md` §2 for every
country result with its interval, and treat any single country figure as weak
evidence on its own.

All four share one 814-feature panel, so they can be scored together in a single
pass. If you retrain one, retrain them all from the same
`prepare_variants.py` run or the tool will refuse to mix them.

## Verifying it works

The standalone predictor reproduces the training pipeline exactly. Scoring 150
isolates from the rifampicin model's own held-out set:

- maximum score difference from the recorded held-out scores: **5e-05**, which is
  the output rounding
- resistant/susceptible calls: **150 of 150 identical**
- against the real laboratory labels for those isolates: accuracy 0.847,
  ROC-AUC 0.947

Run the software tests with:

```bat
.venv\Scripts\python -m unittest discover -s tests
```

69 tests, covering feature-name reproducibility across processes, locus
annotation margins, label handling, metadata normalisation, phenotype merging
and the predictor's matrix construction.

## Retraining or adding a drug

Labels already exist for eleven more drugs in `data/phenotypes.csv`
(streptomycin 1,573; amikacin 1,163; capreomycin 936; and others). To add one:

```bat
.venv\Scripts\python run_full_study.py --drugs STREPTOMYCIN
```

Then add an entry to `DEFAULT_REGISTRY` in `predict_resistance.py`, or pass your
own registry JSON with `--registry`.

## Troubleshooting

| Symptom | Cause |
|---|---|
| "expects different features from" | Models came from different `prepare_variants.py` runs. Retrain them from one panel. |
| Every isolate flagged, all scores low | Positions are not H37Rv coordinates, or the file is empty of resistance-locus variants. |
| "Environment differs from the one these models were fitted in" | Your scikit-learn/pandas/numpy differ from training. Pin `requirements.txt` or retrain. |
| "already exists" | Output paths are never overwritten. Choose a new filename. |
