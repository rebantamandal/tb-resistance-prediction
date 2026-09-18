# Results: TB resistance prediction on the FORUM-TB variant data

Every figure below comes from real laboratory drug susceptibility results joined
to real whole-genome variant calls. Nothing is synthetic, and no label is derived
from the variants used to predict it. All runs are reproducible from §7.

**Research use only.** Not clinically validated, scores uncalibrated, not for
treatment decisions.

---

## The short version

A Random Forest predicting rifampicin resistance from TB genomes scores
**0.958 AUC** under the evaluation protocol most commonly reported, and
**0.686** when tested on a country it was not trained on. The gap is not noise:
the confidence intervals do not overlap. Most of the headline figure is an
artefact of how the data is split, not evidence that the model would work on a
new patient in a new place.

Adding clinical and geographic features — the original planned contribution —
does not reliably help. Under the strictest honest split it makes things
**worse** (−0.039 AUC, interval excludes zero), because those features were
serving as a fingerprint for which study an isolate came from.

---

## 0. The project brief, point by point

| Brief says | Delivered | Status |
|---|---|---|
| Random Forest, for mixed data and interpretable importances | Random Forest, 300 trees; permutation importance used and it caught a real confounder | **Yes** |
| Genomic features: presence/absence of resistance genes | 814 binary variant features across 43 TB resistance loci | **Yes, adapted** |
| Example genes blaTEM, mecA, vanA | Not applicable to M. tuberculosis; the TB equivalents are rpoB, katG, pncA, embB, gyrA and the rest of the 43-locus panel | **Adapted** |
| Clinical: infection site | `infection_site`, from ENA `isolation_source` | **Yes** |
| Clinical: patient age | Not recorded for any of these accessions in any public source | **Not available** |
| Clinical: prior antibiotic exposure | Same | **Not available** |
| Clinical: hospitalisation duration | Same | **Not available** |
| Geographic: country and city | `country` (18 countries) and `city_or_region`, from ENA | **Yes** |
| FORUM-TB genomic-only baseline | Built; the dataset ships no labels, so labels were recovered from BV-BRC, CRyPTIC and NCBI | **Yes, with work** |
| Benchmark ~0.97 rifampicin, ~0.95 isoniazid | 0.958 and 0.947 on independently recovered labels | **Reproduced** |
| Contribution: baseline vs full comparison | Run for four drugs under three evaluation protocols | **Yes** |
| Expected result: full model improves on baseline | It does not. Under study-grouped evaluation it is worse (rifampicin &minus;0.039 AUC, interval excludes zero) | **Negative result** |
| India primary: Delhi, Chennai, Bengaluru | The cohort has 147 Indian isolates, from Tiruvallur (136) and Mumbai (7) | **Not supported** |
| Secondary: USA or UK | USA has 1 isolate; UK has 109 | **UK only, thinly** |
| Third: China or South Africa | South Africa has 1,243 isolates, but 97% are resistant so its held-out metrics are undefined | **Partly** |

Three clinical variables in the brief do not exist for these isolates and were
not invented. The geographic scope had to change because the cohort is what it
is: the countries that can actually carry a comparison are Canada, Malawi, South
Africa and the United Kingdom, not India, the USA and China.

The headline expectation &mdash; that clinical and geographic features would improve on
the genomic baseline &mdash; did not hold. Section 3 explains why, and the explanation
is itself the most useful thing the project found.

---

## 1. Getting labels at all

The FORUM-TB Kaggle dataset publishes one file, `all_variants.csv`, with no
phenotypes — although its description advertises "4 drug resistance labels" and
"247 KB compressed" for what is actually a 555 MB unlabelled variant table
(`DATA_PIPELINE.md` §5). Labels were recovered from three public sources of
*laboratory* DST, joined by sequencing run accession:

| Source | Join route | Contribution |
|---|---|---|
| BV-BRC `genome_amr` | `genome_id` -> registered `sra_accession` | 2,408 isolates |
| CRyPTIC reuse table | `ENA_RUN` | 257 isolates |
| NCBI Pathogen Detection | `Run` -> `AST_phenotypes` | 79 isolates |

Only laboratory evidence was accepted. BV-BRC also publishes its own
classifiers' predictions; training on those would mean learning another model's
output, so they are excluded at the query. Joining BV-BRC through the genome
table's registered accession rather than parsing genome names raised coverage
from 823 to 2,408 isolates.

**2,664 of 9,842 isolates (27%) carry at least one laboratory DST result.**

| Drug | Labelled | Resistant | Susceptible |
|---|---|---|---|
| Rifampicin | 2,625 | 1,446 (55.1%) | 1,179 |
| Isoniazid | 2,603 | 1,994 (76.6%) | 609 |
| Ethambutol | 2,408 | 755 (31.4%) | 1,653 |
| Pyrazinamide | 1,839 | 629 (34.2%) | 1,210 |
| Streptomycin | 1,573 | 914 (58.1%) | 659 |

### The labels themselves check out

Where two sources cover the same isolate and drug, they can be compared:

| Pair | Overlapping results | Agreement |
|---|---|---|
| BV-BRC vs NCBI | 323 | **100%** |
| BV-BRC vs CRyPTIC | 0 | not comparable |
| CRyPTIC vs NCBI | 0 | not comparable |

Perfect agreement on 323 independently sourced results means label noise is not
what limits the numbers below. That matters, because it rules out the easiest
explanation for the generalisation failure in §4.

---

## 2. The evaluation ladder

The same model and the same features, evaluated three ways. Only the split
changes. Rifampicin, 814 binary variant features, 300 trees, threshold 0.50
fixed in advance, 2,000-resample bootstrap for intervals.

| Evaluation | Held out | Baseline AUC | 95% interval |
|---|---|---|---|
| Random split, grouped by **isolate** | 525 | **0.958** | [0.942, 0.973] |
| Random split, grouped by **study** | 701 | **0.908** | [0.885, 0.931] |
| Held-out country: United Kingdom | 93 | **0.758** | [0.620, 0.870] |
| Held-out country: Malawi | 169 | **0.761** | [0.432, 1.000] |
| Held-out country: Canada | 178 | **0.686** | [0.565, 0.799] |

The same collapse happens for every drug tested:

| Drug | By isolate | By study | Held-out countries |
|---|---|---|---|
| Rifampicin | 0.958 | 0.908 | 0.686 &ndash; 0.761 |
| Isoniazid | 0.947 | 0.855 | not run |
| Ethambutol | 0.913 | 0.765 | 0.628 &ndash; 0.843 |
| Pyrazinamide | 0.911 | 0.782 | 0.632 &ndash; 0.837 |

Four drugs, same pattern, no exceptions.

**Row 1 is the number this project set out to reproduce, and it does** — the
brief cited ~0.97 for rifampicin and ~0.95 for isoniazid, and these runs give
0.958 and 0.947 on independently recovered labels. That is a genuine
replication, since the brief's source for those figures was never verifiable.

**Rows 2 to 5 are why that number should not be quoted on its own.**

### Why each rung drops

*Isolate to study (−0.05).* These 9,842 genomes come from **958 different
studies**. Grouping by isolate lets samples from one study sit on both sides of
the split. Each study has its own sampling frame and its own resistance
prevalence, so the model can partly infer the answer from recognising the
cohort. Grouping by study removes that.

*Study to held-out country (−0.2 further).* TB strains cluster into lineages and
outbreaks. Even across studies, a test isolate often has close genetic relatives
in training. Holding out a whole country removes most of that shared population
structure, and performance falls to 0.69–0.76 across three independent
countries. Sensitivity on Canada is 0.559 — barely better than chance on
precisely the cases where a miss harms a patient.

South Africa could not be evaluated: 107 of its 110 labelled isolates are
resistant, so an AUC there would be meaningless. That is reported rather than
quietly dropped.

---

## 3. Baseline versus full model: the planned contribution does not hold

The full model adds real ENA-derived `country`, `city_or_region` and
`infection_site`. Both models are always fitted and scored on **identical rows**.

| Evaluation | Baseline | Full | Difference | 95% interval |
|---|---|---|---|---|
| Random split by study | 0.908 | 0.869 | **−0.039** | [−0.051, −0.028] |
| Held-out Canada | 0.686 | 0.678 | −0.008 | [−0.029, +0.010] |
| Held-out United Kingdom | 0.758 | 0.785 | +0.027 | [+0.002, +0.059] |
| Held-out Malawi | 0.761 | 0.813 | +0.052 | [+0.000, +0.133] |

The sign is not stable. Under the strictest evaluation the extra features make
the model **significantly worse**. On two held-out countries they help slightly,
on one they do not, and no interval is far from zero.

**Conclusion: adding clinical and geographic features has not been shown to
improve prediction, and under study-grouped evaluation it degrades it.** This is
a negative result for the project's original framing, and it is well supported.

### Why: the features were identifying the cohort

An earlier version of the full model included `host` as a clinical feature.
Permutation importance ranked it *second*, above every genomic feature except
rpoB Ser450Leu. That was not biology:

| `host` value | Isolates | Resistant |
|---|---|---|
| `Homo sapiens` | 407 | 13.3% |
| blank | 110 | **97.3%** |

The blank-host isolates were exactly the South African MDR cohort. The model had
learned that *missing metadata* identifies a high-resistance study. `host` was
removed.

`country` does the same thing more subtly. When isolates from one study can
straddle the split, knowing the country helps the model recall that study's
resistance rate, and the full model appears better (+0.017 on the isolate-grouped
run). When study grouping blocks that shortcut, the same features become noise
that costs 0.039 AUC. **The apparent benefit and its disappearance are the same
phenomenon seen from two sides**, which is stronger evidence than either result
alone.

---

## 4. The model is still learning real biology

None of the above means the model is worthless. Permutation importance on the
full model puts `g_rpoB_761155_C_T` — **rpoB Ser450Leu**, the dominant
rifampicin determinant worldwide — far ahead of everything else (AUC drop 0.049,
roughly twice the next feature). katG Ser315Thr and the *fabG1/inhA* promoter
follow. The pipeline recovered the textbook determinants from raw coordinates
without being told what to look for.

The problem is not that the model learned nothing. It is that the causal
mutations alone do not carry performance from 0.96 to a new population, and the
rest of the apparent accuracy was population structure.

---

## 5. What this project should claim

Not: *"we added clinical features and improved on the FORUM-TB benchmark."* The
data does not support it, and the claim would not survive scrutiny.

Instead:

> Published TB resistance benchmarks are typically reported under random splits
> of pooled public genomes. We reproduce such a figure (0.958 AUC for rifampicin)
> and then show that the same model scores 0.69–0.76 on held-out countries, with
> non-overlapping confidence intervals. Grouping the split by study alone costs
> 0.05 AUC. Label noise is excluded as an explanation: two independent sources
> agree on 100% of 323 shared results. We further show that non-genomic metadata
> appears to help under permissive splits and significantly hurts under strict
> ones, because it encodes cohort identity rather than clinical signal.

That is a methodological contribution about evaluation, it is supported by every
number above, and it is more useful than a marginal improvement would have been.

---

## 6. Honest limitations

1. **Three held-out countries, small.** 93–178 isolates each, 7–34 resistant.
   Malawi's interval reaches 1.000. The *pattern* is consistent across all
   three; the individual values are not precise.
2. **Bootstrap resamples rows, not training sets.** These intervals measure
   precision on a given test set. Refitting on different splits would add more
   variance, so treat them as a floor.
3. **27% label coverage**, and only 566 labelled isolates have a country, which
   is what caps the geographic work.
4. **The cohort is MDR-enriched** — katG Ser315Thr is present in 65% of all
   isolates. Prevalence here is far above community TB, so precision and
   sensitivity will not transfer to screening unchanged.
5. **No patient identifiers** anywhere in the source, so repeat isolates from one
   patient may still cross a split even under study grouping.
6. **Variant calls are unfiltered** for depth or quality; the source table
   carries no such fields.
7. **Age, prior antibiotic exposure and hospitalisation duration do not exist**
   in any public record for these accessions and were not invented.

---

## 7. Reproducing everything

```bat
:: metadata and labels
.venv\Scripts\python fetch_sample_metadata.py --accessions data\variant_panel\isolate_variant_matrix.csv --out data\sample_metadata.csv
.venv\Scripts\python fetch_phenotypes.py --out data\phenotype_sources\bvbrc_phenotypes.csv --accessions data\variant_panel\isolate_variant_matrix.csv
.venv\Scripts\python build_phenotype_table.py --bvbrc data\phenotype_sources\bvbrc_phenotypes.csv ^
  --cryptic data\phenotype_sources\CRyPTIC_reuse_table_20240917.csv ^
  --ncbi data\phenotype_sources\ncbi_pathogen_metadata.tsv ^
  --accessions data\variant_panel\isolate_variant_matrix.csv --out data\phenotypes.csv

:: study-grouped matrix, then train
.venv\Scripts\python prepare_variants.py --variants all_variants.csv --out-dir data\rif_bystudy ^
  --min-prevalence 20 --phenotypes data\phenotypes.csv --drug-column RIFAMPICIN ^
  --antibiotic Rifampicin --metadata data\sample_metadata.csv --group-by study
.venv\Scripts\python amr.py train --data data\rif_bystudy\isolate_variant_matrix.csv ^
  --config data\rif_bystudy\config.run.json --out models\rif_bystudy

:: every country held out in turn
.venv\Scripts\python prepare_variants.py --variants all_variants.csv --out-dir data\rif_geo ^
  --min-prevalence 20 --phenotypes data\phenotypes.csv --drug-column RIFAMPICIN ^
  --antibiotic Rifampicin --metadata data\sample_metadata.csv --group-by study --require-geography
.venv\Scripts\python geographic_sweep.py --data data\rif_geo\isolate_variant_matrix.csv ^
  --config data\rif_geo\config.run.json --out-dir models\geo_sweep_rif

:: confidence intervals on any saved run
.venv\Scripts\python analyze_runs.py --runs models\rif_bystudy models\geo_sweep_rif\holdout_Canada
```

New tooling added for this work:

| Script | Purpose |
|---|---|
| `prepare_variants.py` | Long variant table -> isolate x variant matrix; `--group-by study`, `--require-geography` |
| `fetch_sample_metadata.py` | Real country, city and specimen from the ENA portal API |
| `fetch_phenotypes.py` | Laboratory DST from BV-BRC, computational predictions excluded |
| `build_phenotype_table.py` | Merges BV-BRC, CRyPTIC and NCBI with per-label provenance and cross-source agreement |
| `analyze_runs.py` | Paired bootstrap confidence intervals on saved runs |
| `geographic_sweep.py` | One geographic holdout per country, collected into a table |

---

## 8. Next steps, in order of value

1. **Repeat the ladder for isoniazid, ethambutol and pyrazinamide.** Labels
   already exist for all three. If the pattern replicates across four drugs the
   finding is close to unassailable.
2. **Recover more labels.** 27% coverage is the binding constraint, and only 566
   isolates have a country. More labels widen every interval in §2.
3. **Try lineage-aware splitting.** Assigning lineage from the variants and
   holding out whole lineages would separate "population structure" from
   "geography" as explanations, which is currently a reasoned inference rather
   than a measured one.
4. **Report sensitivity at a fixed specificity**, not just AUC. For a clinical
   screening claim that is the more meaningful operating point.
