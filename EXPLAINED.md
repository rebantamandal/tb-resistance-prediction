# The project in plain English

No biology background needed. This explains what the tool does, how to run it,
how to read what comes out, and what the results actually mean.

If you just want the commands, see [`USAGE.md`](USAGE.md). If you want the
numbers and the statistics, see [`RESULTS.md`](RESULTS.md).

---

## 1. What this thing is

**It's a spam filter.**

A spam filter reads an email, looks for suspicious words, and says "spam" or
"not spam."

This reads a bacterium's DNA, looks for suspicious DNA changes, and says "this
antibiotic will fail" or "this antibiotic will work."

That's genuinely the whole idea.

## 2. Why it's worth building

To find out whether a drug works against someone's infection, the standard
method is to grow the bacteria in a lab alongside the drug and see whether it
survives. For tuberculosis that takes **2 to 8 weeks**, because TB grows very
slowly.

Reading the DNA takes **1 to 2 days**.

In the gap between those two numbers, a patient is being treated on a guess. If
the guess is wrong they stay ill and stay infectious. So: can the DNA answer
arrive early enough to be useful? That is the question this project exists to
test.

## 3. What you feed it

One CSV file. It looks like this:

```csv
SAMPLE,CHROM,POS,REF,ALT
Patient1,NC_000962.3,761155,C,T
Patient1,NC_000962.3,2155168,C,G
Patient2,NC_000962.3,1673425,C,T
```

Each row describes **one "typo" in the DNA**. Think of DNA as a very long
document, and this file as a list of spelling errors in it:

| Column | What it means |
|---|---|
| `SAMPLE` | Whose sample this is |
| `CHROM` | Which reference document we're comparing against |
| `POS` | *Where* the typo is (character number 761,155) |
| `REF` | What the letter *should* be (`C`) |
| `ALT` | What it *actually* is (`T`) |

You don't write this file by hand. Whoever runs the DNA sequencing produces it.

### The bit that confuses everyone: rows are not samples

This trips people up, so it's worth being explicit.

**One row is one mutation, not one sample.** A single sample has *thousands* of
mutations, so it takes up thousands of rows.

Here is the example file that ships with this repo:

```
ERR047002     2,496 rows
ERR047009       728 rows
ERR1034590    1,774 rows
ERR1034591    1,783 rows
              ─────
              6,781 rows  ←  but only 4 samples
```

It's shaped like a shopping receipt log:

```
Alice, milk
Alice, bread
Alice, eggs
Bob,   coffee
Bob,   sugar
```

Five rows, two customers. Alice just bought three things. Same idea: 6,781 rows,
four samples, because each sample carries about 1,700 mutations.

If you ran this on 100 patients, you would hand it a file of roughly 170,000
rows and get back 100 answers.

### Why are there so many mutations?

Because **almost none of them have anything to do with drug resistance**. Every
bacterium differs from the reference in thousands of places as ordinary natural
variation — much like any two people's DNA differs in millions of spots without
any of it being medically meaningful.

Out of ~1,700 mutations in a typical sample, only a handful sit in positions
that affect whether a drug works.

**Filtering down to those positions is most of what this project built.** Going
from "here are thousands of mutations" to "here are the 814 positions that
actually matter" is the hard part. The prediction afterwards is comparatively
easy.

## 4. Running it

**The easy way — in a browser.** Run `python predict_app.py`, a page opens, drop
your file on it. No commands to remember, and it highlights the results the tool
is unsure about.

**The command-line way:**

```bat
.venv\Scripts\python predict_resistance.py ^
  --variants examples\sample_isolates.csv ^
  --out test.csv ^
  --trust-local-models
```

You will see something like:

```text
RESEARCH ONLY. Uncalibrated scores, not clinically validated...
Scoring with 5 model(s): Ethambutol, Isoniazid, Pyrazinamide, Rifampicin, Streptomycin
  read 6,781 variant rows, 4 isolates

Wrote test.csv (4 isolates scored).
  Ethambutol         1 resistant /     3 susceptible
  Isoniazid          3 resistant /     1 susceptible
  Pyrazinamide       2 resistant /     2 susceptible
  Rifampicin         2 resistant /     2 susceptible
  Streptomycin       1 resistant /     3 susceptible
```

Line by line:

| Line | Meaning |
|---|---|
| `RESEARCH ONLY...` | Prints every time. Don't treat a real patient with this. |
| `Scoring with 5 model(s)` | Every sample is checked against all five drugs. |
| `read 6,781 variant rows, 4 isolates` | It found 6,781 mutations belonging to 4 samples. |
| `Wrote test.csv` | **The answers are in this file.** |
| `Rifampicin 2 resistant / 2 susceptible` | Just a tally. Not the per-sample answers. |

The terminal only shows a summary. Open the CSV for the real output.

## 5. Reading the output

| Sample | Rifampicin | Isoniazid | Streptomycin | Ethambutol | Pyrazinamide |
|---|---|---|---|---|---|
| ERR047002 | Susceptible `0.25` | **Resistant** `0.74` | Susceptible `0.24` | Susceptible `0.19` | Susceptible `0.16` |
| ERR047009 | Susceptible `0.17` | Susceptible `0.16` | Susceptible `0.27` | Susceptible `0.20` | Susceptible `0.18` |
| ERR1034590 | **Resistant** `0.91` | **Resistant** `0.98` | **Resistant** `0.68` | **Resistant** `0.76` | **Resistant** `0.72` |
| ERR1034591 | **Resistant** `0.85` | **Resistant** `0.89` | Susceptible `0.43` | Susceptible `0.45` | **Resistant** `0.69` |

Read a row as: *"For sample ERR1034590, rifampicin will probably fail."*

- **Resistant** — the drug probably won't work
- **Susceptible** — the drug probably will work
- **The number** — confidence, 0 to 1. At or above 0.50 it says Resistant.

Three things worth noticing in that table:

- **ERR1034590 is bad news across the board.** All five drugs flagged. That is
  what a heavily drug-resistant sample looks like.
- **ERR047009 is clean.** Everything low. Standard drugs should work.
- **ERR1034591's ethambutol is `0.45` and streptomycin `0.43`** — right on the
  fence. Both were called Susceptible only because they sit below 0.50. Anything between about 0.4 and
  0.6 is the tool saying "honestly, not sure." In real use those are the ones
  you would send for laboratory testing rather than trusting.

### The score is not a percentage

`0.90` does **not** mean "90% chance this is resistant." The models were never
calibrated to produce probabilities. Treat the number as a ranking — higher
means more confident — and nothing finer than that.

### The column to always glance at

`known_features_present` tells you how many recognisable mutations were found:

| Sample | Mutations in file | Ones the model uses |
|---|---|---|
| ERR047002 | 2,496 | 41 |
| ERR047009 | 728 | 5 |
| ERR1034590 | 1,774 | 30 |

**If this is `0`, do not read the result as "Susceptible."** It means the tool
recognised nothing, almost always because the file was produced against the
wrong reference. Zero recognised mutations looks identical to a perfectly clean
sample, which is why there is a `warning` column that fills in when it happens.

## 6. Is it any good?

The example above uses four samples whose real laboratory results are known, and
which the model never saw during training. It got **4 out of 4 correct**.

Four samples proves the tool runs, not that it is accurate. For accuracy we
tested it properly — and that produced the most interesting part of the project.

### What "0.958" means

The standard score for this kind of tool is called **AUC**. Ignore the name. It
means:

> Show the tool one resistant sample and one non-resistant sample. How often
> does it correctly say which is which?

- **0.5** = coin flip, useless
- **1.0** = perfect
- **0.958** = right about 96% of the time

The published benchmark this project set out to match was about 0.97. We got
0.958 — matched.

## 7. The surprising part: the score depends on how you grade the exam

Think of the model as a student sitting a test.

You train it on some samples, then test it on samples it has never seen. Sounds
fair. But **which** samples you hold back changes everything.

Here is the catch. Those 9,842 samples came from **958 different collections** —
958 separate hospitals and research projects that each gathered samples for
their own reasons. Each collection has its own character. One was studying a
drug-resistant outbreak, so 97% of its samples are resistant. Another was
routine screening, so almost none are.

If you split the data randomly, samples from the **same hospital** end up in both
the study material and the exam. So the model can cheat — not deliberately, it
just notices:

> "This looks like the samples from Hospital X. Those are nearly always
> resistant. I'll say resistant."

It gets the answer right **without understanding any biology at all.** It is
recognising the hospital, not the disease.

Split it so that no hospital appears on both sides, and the score drops:

| Drug | Split randomly | Split by hospital | Drop |
|---|---|---|---|
| Rifampicin | 0.958 | **0.908** | &minus;0.05 |
| Isoniazid | 0.947 | **0.855** | &minus;0.09 |
| Ethambutol | 0.913 | **0.765** | &minus;0.15 |
| Pyrazinamide | 0.911 | **0.782** | &minus;0.13 |

All four drugs. Every time.

**That gap is how much of the published number comes from recognising the
collection rather than the biology.** The 0.958 is not fake — it just answers an
easier question than the one that matters in practice.

This is the project's strongest finding. **Plan around the middle column.**

## 7b. We checked this on a second, much bigger dataset

A finding measured once could be a quirk of the data. So the whole test was
repeated on a completely different collection of samples — the **CRyPTIC**
dataset, about **seventeen times larger** (45,141 samples instead of 2,625).

**Same result.** Splitting by hospital instead of randomly cost 0.044 for
rifampicin and 0.043 for isoniazid, against 0.050 and 0.092 the first time. The
effect is real, not a fluke of one dataset.

### And it killed our explanation

We had assumed the model recognises hospitals because **samples from one
hospital are genetically related** — same outbreak, same family of bacteria. It
is the obvious guess.

The CRyPTIC data records each sample's genetic family, so we could test it
properly: train the model with **every close relative removed**, then test.

If relatedness were the answer, the score should have dropped. **It didn't
move at all** — 0.954 against 0.952. Same for isoniazid.

So that explanation is wrong, and we only found out by measuring it.

### What is actually going on

The clue is in *which* kind of mistake increases. On an unfamiliar hospital:

- **Catching resistance stays fine** — 0.934 to 0.930, barely changed
- **False alarms jump** — correct "this is fine" calls fall from 0.802 to 0.599

So the model has genuinely learned the resistance mutations. What it also
learned is each hospital's normal background — which *harmless* mutations are
common there — and it uses "none of the usual harmless stuff" as evidence that a
sample is fine. At a new hospital the usual background is different, so healthy
samples start looking suspicious.

**What this means if you use it:** at a new hospital, a "Resistant" result
deserves more doubt than the headline number suggests. A "Susceptible" result is
about as trustworthy as advertised.

### One more check, and it agrees

We also just asked the model directly: *which DNA changes are you actually
using?* The CRyPTIC data labels its mutations with proper names, so the answer
is readable.

Top of the list, by a factor of **seven**, is the single mutation that medicine
has known causes rifampicin resistance for decades. Four of its top eight are
changes in the same small stretch of that one gene.

Meanwhile the changes that merely mark a bacterium's family tree — carried by
about 37,000 of the 45,000 samples — sit near the bottom, twenty times less
important.

So two completely different checks agree: the model is using real biology, not
family resemblance.

One honest wrinkle: two mutations high on its list actually cause resistance to
*other* drugs. That is not a mistake — bacteria resistant to one first-line drug
are usually resistant to several, so those mutations genuinely do predict
rifampicin resistance. But it is guilt by association rather than cause, and it
would mislead the model on an unusual strain resistant to only one drug.

## 8. Something we got wrong, and corrected

We then tried something harder still: remove an entire *country* from training,
then test on it.

Rifampicin dropped a lot, from 0.908 to about 0.69. An earlier version of this
project's write-up concluded "performance collapses across geography."

**Then we checked the other three drugs and it wasn't true.** Isoniazid scored
**0.936** on a country it had never seen — *better* than the 0.855 it managed on
familiar data. Six of seventeen country tests beat their baseline. And streptomycin scored **0.359 on South Africa — worse than a coin flip**, meaning it got that group reliably backwards.

The original claim generalised from a single drug. It is corrected throughout
now. The honest version is: one drug degrades across countries, others do not, and one fails outright on a single group,
and the numbers are too noisy to conclude much either way.

It is recorded here rather than quietly deleted because a reviewer checking the
numbers would have found it, and because "we checked and corrected ourselves" is
a better position than being caught.

## 9. The question that cannot be answered

The project's original plan was: add location and patient information on top of
the DNA, and show the combination beats DNA alone.

We tested it. It gave **+0.017** measured one way and **&minus;0.039** measured
another way. Opposite answers, both apparently significant.

Here is why, and it is clean:

**Every single collection in this dataset happened in exactly one country.** All
the Canadian samples come from one project. All the South African samples come
from one project. No project ever spans a border.

So "which country" and "which hospital" are **the same piece of information**.
There is no way to tell them apart.

> It is like asking whether students do better because of the teacher or because
> of the classroom — when every teacher only ever taught in one classroom. You
> can measure the combination. You can never separate them.

So the answer is not "location doesn't help." It is **"this dataset cannot
answer that question."** That is a legitimate finding, and more defensible than
a number that is secretly measuring something else.

It is also why the shipped models use DNA only.

## 10. What to claim if you present this

**Don't say:** *"We added location data and beat the benchmark."* The data
cannot support it.

**Do say:**

> We reproduced the published benchmark (0.958 AUC for rifampicin). We then
> showed that the standard way of testing these models inflates the score by
> 0.05 to 0.15 across five drugs, because samples from the same collection end
> up on both sides of the test. We also showed that this dataset cannot test
> whether location adds anything, because location and collection are the same
> variable within it.

Every number there is backed by a confidence interval in
[`RESULTS.md`](RESULTS.md), and "here is why the standard benchmark is
optimistic" is a more interesting contribution than "we got +0.02."

## 11. One more point in the project's favour

The model genuinely learned real biology.

Ask it which DNA position mattered most for its rifampicin predictions, and it
points to **the exact mutation that microbiologists have known causes rifampicin
resistance for decades** — roughly twice as important as anything else it uses.

Nobody told it where to look. It found that from raw position numbers.

So it is not only pattern-matching on hospitals. It found the right answer for
the right reason, which is the strongest single piece of evidence that the
pipeline is built correctly.

---

## Where to go next

| Document | What's in it |
|---|---|
| [`USAGE.md`](USAGE.md) | Commands, input format, troubleshooting |
| [`examples/`](examples/) | A file you can run immediately |
| [`RESULTS.md`](RESULTS.md) | Every number, with confidence intervals |
| [`DATA_PIPELINE.md`](DATA_PIPELINE.md) | How the raw data became a usable dataset |
| [`README.md`](README.md) | The overview |
