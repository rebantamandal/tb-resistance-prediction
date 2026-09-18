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
**cannot be evaluated on this cohort at all**. Every study in it is confined to
exactly one country (Cramér's V between study and country = 1.000), so
geography and cohort identity are the same variable. Two earlier comparisons
returned +0.017 and −0.039 AUC; both intervals exclude zero and they disagree in
sign, because each was measuring a different artefact rather than geography. §3
sets this out.

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
| Expected result: full model improves on baseline | Not answerable with this cohort: study and country are perfectly confounded (Cram&eacute;r's V 1.000), and only 5 independent groups have complete metadata | **Untestable** |
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

## 3. Do clinical and geographic features help? The question is not answerable here

This was the project's stated contribution, so it deserves a precise answer
rather than a number. The answer is that **this cohort cannot test the
hypothesis**, and the two earlier attempts each measured a different defect.

### Two comparisons, two opposite answers, both wrong

| Comparison | Full minus baseline | Why it is not informative |
|---|---:|---|
| Complete metadata, grouped by **isolate** | **+0.017** [+0.003, +0.041] | Isolates from one study sit on both sides of the split, so `country` works as a lookup for that study's resistance rate. |
| All labelled isolates, grouped by **study** | **&minus;0.039** [&minus;0.051, &minus;0.028] | `country` is blank for 78% of rows, and blankness itself predicts resistance (61.8% resistant when missing vs 30.6% when present). The model is handed a missingness signal, not a geographic one. |

Both intervals exclude zero, and they have opposite signs. That alone should
prevent either being reported as the project's result.

### Why no third comparison can fix it

The obvious fix is to run the cell neither covered: complete metadata **and**
grouped by study. That cohort has 517 labelled isolates. It also has, in total:

| | Count |
|---|---:|
| Independent studies | **5** |
| Countries | 4 |
| Studies spanning more than one country | **0** |
| Cram&eacute;r's V between study and country | **1.000** |

Across the whole labelled set with a recorded country (570 isolates), **100% of
studies are confined to exactly one country**, and Cram&eacute;r's V is 1.000. Study
identity and country are not merely correlated in this data; they are the same
variable measured twice.

The consequences are unavoidable:

- **Group by isolate** and `country` becomes a cohort lookup, inflating the full
  model. That is the +0.017.
- **Group by study** and every held-out isolate belongs to a country absent from
  training, so `country` is an unseen category carrying no information. The full
  model can only match or trail the baseline, whatever geography's true effect.

There is no split of this dataset under which a genuine geographic effect could
be distinguished from a cohort effect, because no study ever crosses a border.

### The clinical feature has the same problem

`infection_site` varies within only **2 of the 5** studies, both in Malawi. In
the other three it is constant, so there it is a pure study indicator. Within
the two where it does vary, there are **7 resistant isolates in total** &mdash; no
power to detect anything.

### Sample size compounds it

The complete-metadata cohort holds 5 independent groups. A grouped 80/20 split
puts **one study** in the test set. For rifampicin that split is degenerate
(`PRJEB2358` contains no resistant isolate at all), so the run cannot even
complete. Five groups is below any reasonable threshold for a grouped
comparison.

### What can honestly be said

1. The clinical and geographic features available for this cohort **cannot be
   evaluated**, because each is perfectly confounded with study identity.
2. The earlier +0.017 and &minus;0.039 are measurements of that confounding, not of
   geography. Neither belongs in a write-up as a finding about resistance.
3. The released models are therefore **genomic-only**. That is a decision forced
   by the data, not a conclusion that geography is irrelevant to resistance &mdash;
   which is almost certainly false in reality, and remains untested here.

### What would answer the question

A cohort where the two are separable: **isolates from several countries within
the same study**, or **several independent studies within one country**. Either
breaks the confounding. The CRyPTIC compendium collects isolates from 23
countries under one protocol and is the obvious candidate; only 257 of its
isolates overlap this variant file, so it would mean starting from CRyPTIC's own
genomes rather than FORUM-TB's.

Recording the confounding explicitly is worth more than reporting a number that
measures it by accident.

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
