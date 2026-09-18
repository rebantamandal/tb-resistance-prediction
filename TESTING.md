# Verification record

## Checked in this build

- **38 automated tests passed**, using Python's built-in unittest runner. The tests cover grouped splitting, identical holdout rows for both models, training-only imputation statistics, model persistence, single/batch scoring, missing values, invalid labels, invalid genomic codes, negative age, duplicate isolates, out-of-scope organisms/drugs, all-missing inputs, geographic holdouts, undefined subgroup metrics, feature importance, overwrite protection, and local request checks.
- The actual command-line entry points were exercised for training, saving/loading, and batch prediction.
- The browser JavaScript was exercised in Chromium: CSV inspection, column mapping, full/baseline training, run selection, single-sample prediction, three-row batch prediction, CSV download, and a 390-pixel responsive layout. No JavaScript exceptions were observed in that workflow. Desktop views were visually inspected.
- The browser checks used a Python request relay to the real local HTTP server because this container's Chromium policy blocks direct navigation to loopback URLs. The local HTTP endpoints were also tested directly without the browser relay. The relay is a test-environment accommodation, not part of the delivered app.

## Environment actually used

Python 3.13.5; scikit-learn 1.8.0; pandas 2.2.3; NumPy 2.3.5; joblib 1.5.3. The code targets Python 3.11 or newer; the user's Windows machine was not available for direct testing.

## What the tests do NOT establish

All test observations were abstract synthetic fixtures with explicitly non-biological organism/drug names. No real patient dataset was used. No synthetic model, prediction sample, or performance score is included as a project result. Running the tests creates temporary models and then removes them.

These checks do not establish clinical validity, real-world accuracy, geographical transferability, calibration, absence of all data leakage, reproducibility of the claimed FORUM-TB results, or suitability for medical decision-making. Those require real, correctly linked data and an appropriate scientific validation process.

## Re-run

```bat
.venv\Scripts\python -m unittest discover -s tests -v
```

The test fixtures are in `tests/test_model.py`; they are only for software checks, not for training the research project.
