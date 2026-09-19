# Runnable example

`sample_isolates.csv` holds the real variant calls for four *M. tuberculosis*
isolates, taken from the public FORUM-TB data. All four were in the rifampicin
model's **held-out test set**, so the model never trained on them, and their
true laboratory results are known.

## Run it

```bat
.venv\Scripts\python predict_resistance.py ^
  --variants examples\sample_isolates.csv ^
  --out my_predictions.csv ^
  --trust-local-models
```

On macOS or Linux, use `.venv/bin/python` and `\` line continuations.

## What you should see

```text
isolate_id   known_features  Rifampicin   score   Isoniazid    score   Streptomycin  score
ERR1034590   30              Resistant    0.91    Resistant    0.98    Resistant     0.68
ERR1034591   30              Resistant    0.85    Resistant    0.89    Susceptible   0.43
ERR047002    41              Susceptible  0.25    Resistant    0.74    Susceptible   0.24
ERR047009     5              Susceptible  0.17    Susceptible  0.16    Susceptible   0.27
```

All five drugs are scored; only three are shown here for width.

Against the real laboratory rifampicin results, all four calls are correct:

| Isolate | Predicted | Laboratory result |
|---|---|---|
| ERR1034590 | Resistant | Resistant |
| ERR1034591 | Resistant | Resistant |
| ERR047002 | Susceptible | Susceptible |
| ERR047009 | Susceptible | Susceptible |

Four isolates is a demonstration that the tool runs correctly end to end, not a
performance measurement. For that, see [`RESULTS.md`](../RESULTS.md).

`predictions.csv` and `predictions.report.json` in this folder are the committed
output of the run above, so you can diff your result against them.

## Using your own data

Replace `sample_isolates.csv` with your own file in the same five-column format.
Positions must be called against the H37Rv reference (`NC_000962.3`) or the
model will see no features. See [`USAGE.md`](../USAGE.md) for the full spec.
