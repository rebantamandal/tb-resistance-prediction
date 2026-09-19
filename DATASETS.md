# Datasets that could take this project further

This project has three concrete gaps. This document lists the datasets that
would close each one, what is actually in them, and how to get them. Everything
here was checked directly rather than assumed.

## The gaps, in priority order

| # | Gap | Effect on the project |
|---|---|---|
| 1 | **Location and collection are the same variable** | The original research question cannot be answered at all ([`RESULTS.md`](RESULTS.md) §4) |
| 2 | **No patient identifiers** | Repeat samples from one patient can land on both sides of a test, flattering every figure |
| 3 | **No clinical details** | Age, previous treatment and hospital stay don't exist for these samples |

---

## 1. CRyPTIC — the best single option

**What it is:** a research consortium that collected TB samples from 23
countries and tested them all **using one standardised laboratory method**.

**Why it matters here:** the whole reason this project cannot answer its
original question is that every collection in the current data sits in exactly
one country. CRyPTIC deliberately used one protocol across many countries, which
is the structure needed to separate the two.

**Verified contents** (checked against their public FTP):

| File | Size | What's in it |
|---|---|---|
| `GENOMES.csv.gz` | 5 MB | 77,860 genomes, with **patient IDs** (`SUBJID`, 76,492 unique) and **genetic lineage** (`LINEAGE_NAME`) |
| `DST_MEASUREMENTS.csv.gz` | 3.3 MB | Laboratory results, one standardised method |
| `MUTATIONS.csv.gz` | 1.4 GB | Per-sample DNA changes — the equivalent of `all_variants.csv` |
| `MYKROBE_LINEAGE.csv.gz` | 0.5 MB | Lineage assignments |
| `CRyPTIC_reuse_table_*.csv` | 5 MB | 12,287 samples, 13 drugs, ready-made |

**This closes gaps 1 and 2 at once**, and adds something the project currently
lacks entirely: **lineage**.

Lineage matters because the project's main finding is that the model partly
succeeds by recognising which collection a sample came from. The most likely
underlying reason is that samples from one collection are genetically related —
same lineage, same outbreak. With lineage recorded, that becomes testable
directly: hold out a whole lineage and see what happens. Right now it is a
reasoned explanation rather than a measured one.

**Access:** entirely open, no registration.
`https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/`

**Cost:** the 1.4 GB mutations file is the main download. Only 257 of CRyPTIC's
samples overlap the current dataset, so this means working with their genomes
rather than extending the existing ones.

**Recommended first step:** take `GENOMES.csv.gz` (5 MB) and the reuse table and
check how many of the 25 collection sites appear in more than one country. If
even a few do, the geography question becomes answerable.

---

## 2. TB Portals — the only realistic source of clinical data

**What it is:** a US National Institutes of Health programme collecting TB
patient records from nine countries, including **India and South Africa**.

**What's in it:** over 33,000 patient cases with linked clinical, demographic,
laboratory, X-ray and genomic data. One published subset has 2,428 cases from
Azerbaijan, Belarus, Moldova, Georgia, Romania, China, India, Kazakhstan and
South Africa, of which 1,611 (66%) are multi-drug resistant and 952 (39%) have
genomic data.

**Why it matters here:** this is the **only** source found that has the clinical
fields the original project brief asked for — patient age, previous antibiotic
treatment, HIV status, hospital stay — *linked to the same patients' genomes*.
Those fields do not exist in any public record for the current samples, which is
why three rows of the brief are marked "not available".

It also covers India and South Africa directly, which the current data barely
does (India has 147 samples here; the United States has 1).

**Access:** free, but by **request**, not instant download. Genomic and imaging
data through TB Portals; clinical data through NIAID's clinical data access
process. Expect a form and a wait.

- `https://tbportals.niaid.nih.gov/`
- `https://depot.tbportals.niaid.nih.gov` — browse before requesting

**Recommended first step:** browse the DEPOT explorer to confirm the fields are
what the brief needs before committing to the request process.

---

## 3. WHO mutation catalogue — a reference, not training data

**What it is:** the World Health Organization's catalogue of which DNA changes
cause resistance, built from over 52,000 samples across 67 countries for 13
drugs.

**What it's useful for:** *checking* this project's models rather than training
them. The catalogue says which DNA positions are known to matter. Comparing that
against what the models actually rely on is a strong independent validation —
the project already does an informal version of this by noting that the top
feature is the known rifampicin mutation.

It would also let feature names be reported as proper mutation names rather than
position numbers, which would make the output far more readable to a specialist.

**Important limit:** the public release is a catalogue of *mutations*, not a
table of *samples*. It cannot be used to train or test a model — there are no
per-sample records in it.

**Access:** open.
`https://www.who.int/publications/i/item/9789240082410`

---

## 4. Already in use

| Source | Role | Status |
|---|---|---|
| BV-BRC | 2,408 laboratory results | Fully used |
| NCBI Pathogen Detection | 79 results, and a cross-check | Fully used |
| ENA Portal | Country, city, specimen for all 9,842 samples | Fully used |
| CRyPTIC reuse table | 257 results | **Only the phenotypes are used** — the genomes and lineage are not |

Worth noting: the existing pipeline already pulls from CRyPTIC, but only takes
the laboratory results. The genome and lineage tables described in section 1 are
sitting in the same place, unused.

---

## What to do with limited time

1. **Download CRyPTIC's `GENOMES.csv.gz` (5 MB) and check the site-to-country
   mapping.** One afternoon. If sites span countries, the project's original
   question becomes answerable and that changes the whole write-up.
2. **Use CRyPTIC's lineage data to test the population-structure explanation.**
   The project's central claim is that the model recognises collections rather
   than biology. Holding out a whole lineage would test that directly instead of
   inferring it.
3. **Request TB Portals access** if clinical features genuinely matter to the
   project. It is the only route, and the request takes time, so start early.
4. **Map the model's features to WHO mutation names.** Cheap, and it makes the
   output legible to anyone who knows the field.

## What will not help

- **More Kaggle TB datasets.** The one this project uses was checked against
  Kaggle's own listing and is a single unlabelled file. Others found were chest
  X-ray image collections, which are a different problem entirely.
- **Simply finding more samples.** Adding samples from yet more single-country
  collections makes the confounding worse, not better. What is needed is
  *structure* — collections that cross borders — not volume.
