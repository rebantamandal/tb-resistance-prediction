# From `all_variants.csv` to a trainable matrix

This document covers the step that sits between the raw Kaggle download and the
existing `amr.py` trainer. It explains what `prepare_variants.py` does, what it
produced on the supplied data, and the one thing still missing before any model
can be trained.

## 1. What the supplied file actually is

`all_variants.csv` is a long-format variant call table, not a feature matrix:

```text
SAMPLE,CHROM,POS,REF,ALT
ERR038266,NC_000962.3,1849,C,A
```

| Property | Value |
|---|---|
| Rows | 15,057,917 variant calls |
| Isolates (`SAMPLE`) | 9,842 ENA/SRA run accessions |
| Reference contig | `NC_000962.3` only (M. tuberculosis H37Rv) |
| Distinct variants | 758,768 |
| Median calls per isolate | ~1,530 |
| Resistance phenotype | **absent** |

The variant distribution is dominated by rare calls: 567,188 variants (75%)
occur in exactly one isolate, and only 4,443 occur in 500 or more. Fitting a
Random Forest to the untrimmed 758k-column space would mostly be fitting
sequencing noise and lineage background, so the preparation step filters.

`Downloads\archive\all_variants.csv` is a byte-identical copy of the same file;
only one is needed.

## 2. What `prepare_variants.py` does

Two streaming passes over the file, so peak memory stays independent of its
555 MB size.

**Pass 1 — prevalence and annotation.** Every variant is keyed as
`<pos>_<ref>_<alt>` and counted by the number of *distinct isolates* carrying
it. Each position is mapped to an H37Rv gene using
`resources/tb_resistance_loci.csv`, a curated table of 43 loci with
drug associations. Each window is widened by `--margin` bases (default 200) on
both sides, because several decisive determinants are promoter variants outside
the coding sequence — most importantly the *fabG1/inhA* promoter.

**Pass 2 — matrix fill.** The retained variants become binary columns and the
matrix is filled isolate by isolate.

Two filters control dimensionality:

- `--min-prevalence` (default 20 isolates) drops the singleton tail.
- `--max-prevalence-fraction` (default 0.99) drops near-fixed sites, which
  carry essentially no discriminative signal.

`--panel resistance` (default) keeps only variants inside the locus table, which
is what makes the resulting feature importances readable as gene-level
statements. `--panel genome-wide` keeps everything passing the prevalence
filter, for a comparison run that does not presuppose the candidate gene list.

## 3. What it produced

```bat
.venv\Scripts\python prepare_variants.py ^
  --variants all_variants.csv --out-dir data\variant_panel --min-prevalence 20
```

**9,842 isolates x 814 variant features** across 42 genes, from 17,990 distinct
variants found inside the locus windows. Output in `data\variant_panel\`:

| File | Content |
|---|---|
| `isolate_variant_matrix.csv` | The matrix, in the column layout `amr.py` expects (16.5 MB). |
| `feature_dictionary.csv` | Every feature with its position, alleles, gene, associated drugs and carrier count. |
| `config.generated.json` | A ready-to-run `amr.py` config listing all 814 genomic columns. |
| `prep_report.json` | Parameters, counts, per-gene feature totals, and limitations. |

Genes contributing the most features: `rrs` (93), `rrl` (88), `pncA` (55),
`rpoC` (53), `rpoB` (50), `embB` (46), `gid` (42), `gyrA` (34).

The matrix is 16.5 MB, slightly above the browser uploader's 15 MiB cap, so
train it through the command line. Raising `--min-prevalence` brings it under
the cap if the browser workflow is preferred.

### The output is biologically coherent

The highest-prevalence features are the textbook resistance determinants, which
is the strongest available evidence that coordinates, alleles and the annotation
margin are all correct:

| Feature | Isolates | Identity |
|---|---|---|
| `g_katG_2155168_C_G` | 6,404 (65%) | katG Ser315Thr, the dominant isoniazid determinant |
| `g_rpoB_761155_C_T` | 4,425 (45%) | rpoB Ser450Leu, the dominant rifampicin determinant |
| `g_fabG1_1673425_C_T` | 1,565 | *fabG1/inhA* promoter -15C>T (isoniazid, ethionamide) |
| `g_rpoB_761139_C_T` / `_C_G` | 306 / 285 | RRDR codon 445 hotspot |

Those frequencies also characterise the cohort: it is heavily enriched for
drug-resistant TB rather than being a population sample. Resistance prevalence
in any model trained on it will be far above community prevalence, so precision
and sensitivity will not transfer to a screening setting unchanged.

## 4. Real clinical and geographic features from ENA

The isolate identifiers are ENA/SRA run accessions, so the samples have publicly
registered metadata. `fetch_sample_metadata.py` retrieves it from the ENA portal
API, which is public and needs no account:

```bat
.venv\Scripts\python fetch_sample_metadata.py ^
  --accessions data\variant_panel\isolate_variant_matrix.csv --out data\sample_metadata.csv
```

All 9,842 accessions were found. Coverage is partial, because it is only ever
whatever the original submitter recorded:

| Field | Isolates | Feature category |
|---|---|---|
| `country` | 4,628 across 18 countries | geographic |
| `city_or_region` | 2,175 | geographic |
| `infection_site` (from `isolation_source`) | 3,682 | clinical |
| `collection_year` | 4,796 (1995-2016) | temporal |
| Country **and** infection site | **3,250** | both |

Countries: South Africa 1,243, Peru 889, Russia 593, Belarus 468, Sweden 415,
Canada 179, Malawi 170, Mali 164, **India 147**, Thailand 126, United Kingdom
109, Uganda 47, Romania 34, and five smaller. Cities include Lima (889),
Bangkok (126), **Tiruvallur (136)** and Mumbai (7) in India, Durban, Oxford,
Kampala and British Columbia.

Specimen categories: respiratory 1,309, unspecified 1,199, unspecified fluid
1,061, extrapulmonary 73, other 40.

### Normalisation is explicit, not guessed

ENA writes the literal string `missing` into these fields, and some submitters
put a hospital name in the country field. `resources/metadata_normalization.json`
holds the rules, so every decision is auditable and editable:

- Null placeholders (`missing`, `not known`, `Not applicable`, ...) become blank.
  Before this step the country field showed 5,799 "populated" values; 1,091 of
  them were the word `missing`.
- 64 isolates name a KwaZulu-Natal hospital instead of a country. These are
  **not** recoded to South Africa, because that would be an inference rather than
  a recorded fact. The verbatim value is kept in `collection_site` and the
  country left blank. Recode them only after verifying against the source study.
- `isolation_source` free text is grouped into coarse categories. Strings naming
  no anatomical site (`patient`, `culture`, `clinical sample`) map to
  `unspecified`, which is a recorded non-answer and must not be read as
  pulmonary.

### Building the enriched matrix

```bat
.venv\Scripts\python prepare_variants.py --variants all_variants.csv ^
  --out-dir data\panel_with_metadata --min-prevalence 20 ^
  --metadata data\sample_metadata.csv --require-metadata
```

This produced **3,250 isolates x 814 variant features** plus `infection_site`,
`host`, `country` and `city_or_region`. `--require-metadata` keeps only isolates
having both a country and an infection site, so the baseline and full models are
compared on identical rows rather than on different cohorts. The generated
config sets `compare_full: true` and was verified to pass `amr.py`'s own
validation, giving 814 genomic features for the baseline and 814 + 4 for the
full model.

`collection_year` is written into the matrix but deliberately left out of the
generated config. It is a temporal covariate, not a clinical one, and quietly
filing it under "clinical" would misdescribe the feature blocks. Add it to
`clinical_numeric_columns` by hand if you want it, and say so in the write-up.

### This changes the geographic plan in the brief

The brief proposes India as primary (Delhi, Chennai, Bengaluru), USA or UK as
secondary, and China or South Africa as an optional third. This cohort does not
support that:

- **India has 147 isolates**, and the city is Tiruvallur (136), not Delhi,
  Chennai or Bengaluru. Tiruvallur is a genuine Tamil Nadu study site, so it is
  usable, but it is one district cohort rather than three surveillance cities.
- **USA has 1 isolate**; the UK has 109. Neither supports a validation arm.
- The countries that *can* carry a geographic comparison are South Africa
  (1,243), Peru (889), Russia (593) and Belarus (468).

Either re-scope the geographic comparison to South Africa / Peru / Russia, which
the data supports, or keep India in the study and state plainly that it is a
147-isolate single-site subgroup. The `--split-method geography_holdout` option
in `amr.py` works either way. What will not work is presenting an ICMR
three-city Indian analysis on the basis of this file.

## 5. The blocker: no phenotype labels

The Kaggle dataset `nanzhen/forum-tb` contains `all_variants.csv` and nothing
else. There is no per-isolate drug susceptibility testing result, so **no
supervised model can be trained yet**. The pipeline stops at an unlabelled
matrix by design.

### The uploaded file does not match its own description

This was confirmed against Kaggle's public API, which needs no login. The
dataset page describes:

> 9,798 isolates - 2,693 genomic features (AMR gene positions) - 4 drug
> resistance labels (Rifampicin, Isoniazid, Ethambutol, and Pyrazinamide) -
> **247 KB compressed**

The file that is actually published is 555 MB, holds 9,842 isolates, carries
758,768 genome-wide variants rather than 2,693 positions across nine genes, and
contains no labels at all. The file listing confirms `all_variants.csv` is the
only file in the dataset.

The most likely explanation is that the author uploaded an intermediate raw
variant call table instead of the finished ML-ready matrix. The described
version is the one the project brief was planning around, and it does not exist
at that URL. Asking the author on the dataset's discussion tab to publish the
labelled matrix is the most direct route to the phenotypes; the accession list
here is the fallback.

The one shortcut that must not be taken is deriving labels from the variants —
calling an isolate rifampicin-resistant because it carries a known *rpoB*
mutation, then training a model on *rpoB* mutations. That is circular: it would
score near-perfect AUC while measuring only the agreement of a rule with itself,
and it would reproduce none of the real phenotypic discordance that makes this
problem worth modelling.

What is needed is a table of one row per isolate with a documented DST outcome:

```text
SAMPLE,RIFAMPICIN,ISONIAZID
ERR038266,S,S
ERR038737,R,R
```

Because the identifiers are ENA/SRA run accessions, phenotypes can be recovered
from the public sources those accessions came from — the CRyPTIC consortium
release, the Walker et al. and ReSeqTB collections, or the accompanying file on
the Kaggle dataset page if one exists behind the login.

Once that file is in hand:

```bat
.venv\Scripts\python prepare_variants.py ^
  --variants all_variants.csv --out-dir data\rifampicin ^
  --phenotypes phenotypes.csv --drug-column RIFAMPICIN --antibiotic Rifampicin

.venv\Scripts\python amr.py train --data data\rifampicin\isolate_variant_matrix.csv ^
  --config data\rifampicin\config.generated.json --out models\rif_baseline
```

Label handling follows the policy `amr.py` already enforces: `R`/`resistant`/`1`
become the positive class, `S`/`susceptible`/`sensitive`/`0` the negative, and
intermediate, unknown, blank or any other value is **dropped and counted**, not
recoded as susceptible. Conflicting duplicate results for one isolate abort the
run rather than being averaged.

## 6. Verification performed

- The pipeline was run end to end on a 255-isolate subset with a deliberately
  randomised placeholder label column, through `prepare_variants.py` and then
  `amr.py train`. It completed and produced ROC-AUC 0.504. That number is a
  **negative control, not a result**: labels independent of the genotype should
  score at chance, and a pipeline leaking the target would not. It is evidence
  about the software, not about resistance prediction.
- Intermediate (`I`) labels in that fixture were dropped, not recoded: 255
  isolates in, 231 retained, 24 counted as unusable in `prep_report.json`.
- The full two-pass run over all 15,057,917 rows completed and the retained
  features were checked against known H37Rv resistance coordinates, as above.

No performance figure for resistance prediction exists yet, and none can until
real phenotypes are joined.

## 7. Known limitations carried into any model

These are also written into `prep_report.json` for every run.

1. **No patient identifier.** The source has none, so `patient_id` repeats
   `isolate_id`. If the cohort contains multiple isolates per patient, related
   samples can land on both sides of the train/test split and inflate measured
   performance. This is the single most likely source of optimism in the
   baseline and should be stated in any write-up.
2. **No call quality filtering.** The table carries no depth, genotype quality
   or allele fraction fields, so every call is taken at face value.
3. **Coordinate-level, not protein-level.** A variant inside *rpoB* is not
   thereby a resistance-conferring mutation; synonymous and lineage-marker
   variants are retained alongside real determinants. Mapping features to amino
   acid changes against the WHO mutation catalogue would sharpen interpretation.
4. **Cohort enrichment.** As noted above, this is not a population sample.
5. **Partial non-genomic coverage.** Country and infection site exist for
   3,250 of 9,842 isolates, so the full-versus-baseline comparison runs on a
   third of the cohort. Age, prior antibiotic exposure and hospitalisation
   duration are not registered in ENA at all and remain unavailable; they
   cannot be invented.
6. **Metadata is submitter-reported.** Country is where the sample was
   collected, not the patient's residence or country of infection, and the
   free-text fields were never harmonised at source.
