# Results

Every number here comes from real laboratory test results matched to real DNA
data. Nothing is simulated.

Written to be readable without a biology or statistics background. Terms are
explained where they first appear. For the plain-English version of the whole
project, see [`EXPLAINED.md`](EXPLAINED.md).

> **Research use only.** Not medically approved. Scores are not probabilities.
> Not for treatment decisions.

---

## The one-paragraph summary

The tool predicts, from a bacterium's DNA, whether an antibiotic will work.
Tested the way these tools are normally tested, it is right about **95% of the
time**. Tested more strictly — on samples from hospitals it never learned from —
it is right about **85 to 90% of the time**. That gap holds for all five drugs,
and it is the project's main finding: **the standard way of measuring these
tools flatters them.**

That finding has since been **confirmed on a second, independent dataset
seventeen times larger** (section 4), which also rules out the most obvious
explanation for it.

A second question the project set out to answer — whether adding a patient's
location improves predictions — turned out to be **unanswerable with this
dataset**, for a reason explained in section 5.

---

## 1. What the score means

The standard measure is called **AUC**. Ignore the name:

> Show the tool one resistant sample and one non-resistant sample. How often
> does it correctly say which is which?

- **0.5** = coin flip, no better than guessing
- **1.0** = perfect
- **0.95** = right about 95% of the time

Anything **below 0.5 is worse than guessing** — the tool is systematically
getting it backwards. That happens once in this project, in section 3, and it is
reported rather than hidden.

Numbers in `[brackets]` are the **margin of error**. `0.908 [0.885, 0.931]`
means the true value is probably somewhere between 0.885 and 0.931. When two
ranges don't overlap, the difference between them is real rather than luck.

---

## 2. Where the answers came from

The dataset this project was given contained DNA data for 9,842 samples and
**no answers**. Its download page advertised answers for four drugs; the file
actually published had none.

So the answers were recovered from three public medical databases by matching
sample ID numbers:

| Source | Samples matched |
|---|---:|
| BV-BRC (laboratory results only) | 2,408 |
| CRyPTIC research consortium | 257 |
| NCBI Pathogen Detection | 79 |

**2,664 of 9,842 samples (27%) ended up with a real laboratory result.**

Two shortcuts were deliberately avoided:

- **Computer predictions were excluded.** BV-BRC also stores guesses made by
  other software. Training on those would mean learning another program's
  opinions rather than reality.
- **Answers were never derived from the DNA.** The tempting shortcut is to label
  a sample "resistant" because it carries a known resistance mutation, then
  train a model to spot that mutation. That is circular — it would score near
  100% while proving nothing.

### The answers were checked against each other

Where two independent databases covered the same sample and drug, they can be
compared. **They agreed on 323 out of 323.**

That matters, because it rules out "the answers were wrong" as an explanation
for anything below.

### How many samples per drug

| Drug | Samples with an answer | Resistant | Not resistant |
|---|---:|---:|---:|
| Rifampicin | 2,625 | 1,446 | 1,179 |
| Isoniazid | 2,603 | 1,994 | 609 |
| Ethambutol | 2,408 | 755 | 1,653 |
| Pyrazinamide | 1,839 | 629 | 1,210 |
| Streptomycin | 1,573 | 914 | 659 |

---

## 3. How accurate is it?

The same models, tested three ways. **Only the testing method changes** — the
models themselves are identical throughout.

### Test 1 vs Test 2: the main finding

| Drug | Familiar hospitals | **New hospital** | Difference |
|---|---|---|---|
| Rifampicin | 0.958 [0.942, 0.973] | **0.908** [0.885, 0.931] | −0.050 |
| Streptomycin | 0.948 [0.925, 0.969] | **0.897** [0.878, 0.915] | −0.051 |
| Isoniazid | 0.947 [0.928, 0.964] | **0.855** [0.822, 0.887] | −0.092 |
| Ethambutol | 0.913 [0.885, 0.939] | **0.765** [0.723, 0.803] | −0.148 |
| Pyrazinamide | 0.911 [0.877, 0.940] | **0.782** [0.746, 0.819] | −0.129 |

**Five drugs. Every one drops. No margin of error overlaps.**

Here is why it happens. The 9,842 samples came from **958 separate collections**
— different hospitals and research projects that each gathered samples for their
own reasons. Each has its own character: one was studying a drug-resistant
outbreak, so 97% of its samples are resistant; another was routine screening, so
almost none are.

Test 1 splits the samples randomly, so samples from the **same hospital** land
in both the training material and the test. The model can then succeed by
noticing "this looks like the samples from Hospital X, and those are nearly
always resistant" — getting the right answer without understanding anything
about the disease.

Test 2 makes sure no hospital appears on both sides, removing that shortcut.

**The gap is how much of the published figure comes from recognising the
hospital rather than the biology.** Column 1 is the figure normally reported in
the literature. **Column 2 is what to expect in real use.**

The project brief cited roughly 0.97 for rifampicin and 0.95 for isoniazid as
the targets to match. Column 1 gives 0.958 and 0.947 — matched. Column 2 is the
more honest reading of the same models.

### Test 3: removing a whole country

| Drug | New hospital | Results per country held out |
|---|---|---|
| Rifampicin | 0.908 | Canada 0.686 · UK 0.758 · Malawi 0.761 |
| Isoniazid | 0.855 | **Canada 0.936** · **Malawi 0.904** · S. Africa 0.795 · UK 0.773 |
| Streptomycin | 0.897 | Canada 0.856 · Malawi 0.761 · **S. Africa 0.359** |
| Ethambutol | 0.765 | **UK 0.843** · Canada 0.739 · S. Africa 0.628 |
| Pyrazinamide | 0.782 | **Canada 0.837** · S. Africa 0.632 · UK 0.644 |

Bold entries scored **better** than their own column-2 figure.

This is **not** a clean story, and an earlier draft of this document claimed it
was. The honest reading:

- **Rifampicin genuinely degrades** on every country tested.
- **Isoniazid does the opposite** — 0.936 on a country it had never seen, better
  than the 0.855 it manages on familiar data.
- **Streptomycin on South Africa scored 0.359, worse than guessing.** Its margin
  of error is [0.227, 0.493], entirely below 0.5, so this is not a fluke: on
  that group the model is reliably backwards. Those samples are 70% resistant
  and come from a single drug-resistant-outbreak study, so the patterns learned
  elsewhere actively mislead it there.
- **The margins of error are large.** Several span 0.3 or more. Most single
  country results are too uncertain to conclude much from on their own.

**Six of seventeen country results beat their own baseline.** The earlier claim
that "performance collapses across geography" generalised from one drug and was
wrong. It is corrected here rather than quietly deleted, because anyone checking
the numbers would find it.

So: **test 1 to test 2 is a real, consistent, measurable gap.** Test 3 is noisy,
and the streptomycin result is a genuine warning that the tool can fail badly on
an unfamiliar population.

---

## 4. Confirmed on a second, much larger dataset

Section 3's finding was measured on one dataset. To check it was not a quirk of
that data, the whole thing was repeated on the **CRyPTIC consortium** collection
— a completely separate set of samples, roughly **seventeen times larger**, and
collected under one shared laboratory protocol rather than 958 different ones.

| | This project's data | CRyPTIC |
|---|---:|---:|
| Rifampicin samples with an answer | 2,625 | **45,141** |
| Isoniazid samples with an answer | 2,603 | **44,880** |
| Patient identifiers | none | **43,361** |
| Genetic lineage | not recorded | recorded for 71% |

### The main finding replicates

| Drug | Grouped by patient | **Grouped by collection** | Cost |
|---|---|---|---|
| Rifampicin | 0.952 [0.947, 0.956] | **0.908** [0.901, 0.916] | &minus;0.044 |
| Isoniazid | 0.957 [0.953, 0.961] | **0.914** [0.908, 0.921] | &minus;0.043 |

The original measurement on this project's own data was &minus;0.050 for
rifampicin and &minus;0.092 for isoniazid. The direction, the rough size and the
non-overlapping margins of error all hold on independent data seventeen times
the size. **This is no longer a property of one dataset.**

Note also that CRyPTIC allows a *patient*-grouped split, which this project's
data cannot support at all. Even with each patient's samples kept together —
removing one source of flattery entirely — moving to collection-grouping still
costs the same ~0.04.

### And it rules out the obvious explanation

The natural guess is that samples from one collection are **genetically
related** — same outbreak, same family of strain — so a model can recognise the
family rather than the disease. CRyPTIC records each sample's genetic lineage,
which makes that testable rather than assumed.

Holding out **entire genetic lineages**:

| Drug | Grouped by patient | Grouped by **lineage** | Cost |
|---|---|---|---|
| Rifampicin | 0.952 [0.947, 0.956] | 0.954 [0.950, 0.957] | **none** |
| Isoniazid | 0.957 [0.953, 0.961] | 0.959 [0.955, 0.962] | **none** |

**No penalty at all.** The model can be denied every close genetic relative of
the samples it is tested on and lose nothing. Whatever the collection effect is,
**it is not genetic relatedness.**

An earlier draft of this project's write-up offered exactly that explanation.
It is wrong, and it took having lineage data to find that out.

### What the collection effect actually looks like

The sensitivity and specificity split tells the story:

| Drug | Grouping | Sensitivity | Specificity |
|---|---|---:|---:|
| Rifampicin | by patient | 0.934 | 0.802 |
| Rifampicin | **by collection** | 0.930 | **0.599** |
| Isoniazid | by patient | 0.944 | 0.820 |
| Isoniazid | **by collection** | 0.937 | **0.659** |

**Sensitivity barely moves. Specificity collapses.**

On a collection it has never seen, the model still finds nearly every resistant
sample — it has genuinely learned the resistance mutations. What breaks is false
alarms: it flags far more susceptible samples as resistant.

That is the signature of a model that has learned each collection's background —
which harmless mutations are common there — and uses their absence as evidence
of susceptibility. At a new collection those background cues are wrong, so
susceptible samples start looking suspicious.

**Practical consequence:** on a new site, trust a "Resistant" call less than the
headline accuracy suggests, and trust a "Susceptible" call about as much.

### A second, independent check agrees

The lineage test says the model does not lean on genetic family. Asking the
model directly which mutations it depends on says the same thing, by a different
route.

CRyPTIC names its mutations, so the answer is readable:

| Mutation | How much the model depends on it |
|---|---:|
| **rpoB@S450L** | **0.0732** |
| embB@M306V | 0.0108 |
| katG@S315T | 0.0103 |
| embB@M306I | 0.0066 |
| rpoB@H445Y | 0.0059 |
| rpoB@D435V | 0.0054 |
| whiB6@-74_indel | 0.0042 |
| rpoB@H445D | 0.0040 |
| gyrA@E21Q | 0.0034 |
| mmpL5@I948V | 0.0033 |

Three things to read here:

- **rpoB S450L dominates at roughly seven times the next feature.** That is the
  single mutation medicine has known causes rifampicin resistance for decades,
  and the model found it unaided. Four of the top eight are rpoB mutations, all
  clustered in the same small region of that gene.
- **katG S315T and embB M306V rank high but are isoniazid and ethambutol
  mutations**, not rifampicin ones. That is not confusion: strains resistant to
  one first-line drug are often resistant to several, so carrying the isoniazid
  mutation genuinely predicts rifampicin resistance. It is a real correlation,
  but it is co-resistance rather than cause, and a model leaning on it would
  mispredict an unusual strain resistant to only one drug.
- **The lineage markers score near the bottom.** `gyrA@E21Q` and
  `mmpL5@I948V` are carried by roughly 37,000 of the 45,141 samples and mark
  genetic family rather than resistance. The model depends on them about twenty
  times less than on rpoB S450L. Two separate methods — holding lineages out,
  and asking the model what it uses — agree that phylogeny is not what is
  driving the score.

---

## 5. The question that cannot be answered

The project's original plan was to add each patient's **country and clinical
details** on top of the DNA and show the combination does better.

It was measured twice and gave **opposite answers**: +0.017 one way, −0.039 the
other, both apparently meaningful.

Here is why, and it is clean:

**Every collection in this dataset took place in exactly one country.** All the
Canadian samples come from a single project. All the South African samples come
from a single project. No project ever crosses a border.

So **"which country" and "which hospital" are the same piece of information.**
There is no way to separate them.

> It is like asking whether students do better because of the teacher or because
> of the classroom — when every teacher only ever taught in one classroom. You
> can measure the combination. You can never split it apart.

Statistically, the link between "collection" and "country" measures **1.000 on a
0-to-1 scale** — perfectly locked together.

The two earlier measurements were therefore measuring different mistakes:

| Measurement | Result | What it actually captured |
|---|---|---|
| Split randomly | +0.017 | Country acting as a lookup for a hospital's resistance rate |
| Split by hospital | −0.039 | Country left blank for 78% of samples, and *being blank* itself predicts resistance |

Attempting the clean version — complete information *and* split by hospital —
fails for a third reason: only **5 collections** have complete information, and
3 of the 5 drugs cannot be tested at all because the single held-out collection
contains only one type of answer.

**The conclusion is "this dataset cannot answer the question," not "location
doesn't help."** That is a legitimate finding, and more defensible than a number
secretly measuring something else. It is also why the shipped models use DNA
only.

### What would answer it

A dataset where the two come apart: samples from several countries within one
study, or several independent studies within one country. The CRyPTIC consortium
collected samples from 23 countries under a single protocol and would work — but
only 257 of its samples overlap this dataset, so it would mean starting over
with their data instead.

---

## 6. The tool did learn real biology

This is the strongest evidence the pipeline is built correctly.

Ask the rifampicin model which DNA position mattered most to its decisions, and
it points to **the exact mutation that medicine has known causes rifampicin
resistance for decades** — roughly twice as important as anything else it uses.
The next most important are likewise the known isoniazid mutations.

Nobody told it where to look. It found them from raw position numbers.

On the CRyPTIC data the same check is far sharper, because mutations there carry
their proper names: rpoB S450L comes out **seven times** more important than
anything else, and four of the top eight are mutations in the same small region
of rpoB. Section 4 has the full table and the caveat about co-resistance.

So it is not purely recognising hospitals. It found the right answer for the
right reason.

---

## 7. What to claim

**Don't claim:** *"We added location data and beat the benchmark."* The data
cannot support it.

**Don't claim:** *"Performance collapses across geography."* True for one drug,
contradicted by three others.

**Do claim:**

> We reproduced the published benchmark (0.958 for rifampicin, 0.947 for
> isoniazid) using laboratory answers we recovered ourselves. We then showed
> that the standard testing method overstates accuracy by 0.05 to 0.15 across
> five drugs, because samples from the same collection appear on both sides of
> the test. We confirmed this on the CRyPTIC collection — independent data,
> seventeen times larger — where the same effect measures 0.043 to 0.044 for
> both drugs tested. Using CRyPTIC's lineage records we ruled out genetic
> relatedness as the cause: holding out entire lineages costs nothing. The loss
> falls almost entirely on specificity rather than sensitivity, meaning the
> model still finds resistance on an unfamiliar collection but raises far more
> false alarms. Faulty answers are ruled out: two independent sources agree on
> 100% of 323 shared results. We also showed this dataset cannot test whether
> location helps, because location and collection are the same variable within
> it.

Every clause there has a margin of error behind it.

---

## 8. Honest limitations

1. **Scores are not probabilities.** 0.80 does not mean an 80% chance. Use them
   as rankings.
2. **The samples are unusual.** Roughly two thirds carry a major resistance
   mutation. Real-world populations are far less resistant, so accuracy on
   ordinary patients is untested.
3. **It can fail badly on an unfamiliar population** — see streptomycin on South
   Africa scoring below chance.
4. **Only 27% of samples had answers**, which limits how precise every figure
   above can be.
5. **Some margins of error are very wide**, particularly in the country tests.
6. **No patient IDs exist** in the data, so repeat samples from one patient may
   appear on both sides of a test, flattering every figure slightly.
7. **DNA quality is not checked.** The source file carries no quality
   information, so every reported DNA change is taken at face value.
8. **Most patient details do not exist.** Age, previous antibiotic treatment and
   hospital stay length are not recorded anywhere public for these samples, and
   were not invented.

---

## 9. Reproducing this

```bat
:: 1. turn the raw DNA file into a table the model can use
.venv\Scripts\python prepare_variants.py --variants all_variants.csv ^
  --out-dir data\variant_panel --min-prevalence 20

:: 2. collect sample details and laboratory answers from public databases
.venv\Scripts\python fetch_sample_metadata.py ^
  --accessions data\variant_panel\isolate_variant_matrix.csv --out data\sample_metadata.csv
.venv\Scripts\python fetch_phenotypes.py --out data\phenotype_sources\bvbrc_phenotypes.csv ^
  --accessions data\variant_panel\isolate_variant_matrix.csv
.venv\Scripts\python build_phenotype_table.py ^
  --bvbrc data\phenotype_sources\bvbrc_phenotypes.csv ^
  --accessions data\variant_panel\isolate_variant_matrix.csv --out data\phenotypes.csv

:: 3. train and test one drug, all three ways
.venv\Scripts\python run_full_study.py --drugs RIFAMPICIN

:: 4. margins of error for any saved result
.venv\Scripts\python analyze_runs.py --runs models\rif_genomic_baseline
```

Working folders under `data/` are intentionally not stored in this repository.
The commands above rebuild them in a few minutes each.

---

## 10. What to do next

1. **Get more answers.** 27% coverage is the single biggest constraint; more
   would narrow every margin of error here.
2. **Investigate the streptomycin / South Africa failure.** A result below
   chance is unusual and worth understanding before anyone relies on the tool.
3. **Turn scores into real probabilities**, so a number can be read as a
   percentage rather than a ranking.
4. **Report accuracy at a fixed error rate**, not just AUC, which is what
   matters if the tool were ever used as a screening test.
