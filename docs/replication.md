# Replicating the 2014 package

Search has two readings of the SEARCH R package (Haider and Buffa, version 0.0.2, 2014). The
`search2014` preset of the engine computes the method as it was meant, without the defects of
the package ([Porting from R](porting.md)). This page is about the other reading:
`search.r2014`, a port that reproduces the package as it was written, defects included, and
about the question it was built to answer: **can a result of the 2014 package be reproduced
exactly?**

The study was run on 2026-11-19 with Python 3.12.9, NumPy 2.5.3 and SciPy 1.18.1 on macOS
(arm64), against two installations of R: R 4.4.3 on the same machine, and R 3.1.0 for Linux
x86-64, the R released three days after the date of the package, run from the container image
`rocker/r-ver:3.1.0`. On 2026-11-20 every function of the package was compared with R on its
own; that comparison found seven differences the whole runs had not reached, and they were
corrected ([Function by function](#function-by-function)).

## The answer

- **Yes, for everything the package does with correlations**, which is its default and its
  documented use. `search.r2014` returns what the unmodified package returns in R 4.4.3: the
  same genes under the same names in the same order, every number the same double bit for
  bit, the same files character for character, and where R stops with an error, the same
  error at the same point. This holds in all 22,535 cases compared (37.9 million numbers),
  among them the examples of the manual, every known defect, bootstrap and jackknife runs
  with their summary, and three growth iterations on the GBM IT cohort (21254 genes).
- **Function by function as well.** Each of the 14 ported functions of the package, called on
  its own in R 4.4.3 and in the port, returns the same in every case tried, and so does each
  of the 18 rules of R the port had to write down. Getting there took seven corrections to
  the port, all in states no whole run had reached.
- **Yes also for the R of 2014, as far as it was run.** In R 3.1.0 on x86-64 the 26 cases on
  the package's own data are reproduced in the same exact sense, once the port is told to
  compute as that R does.
- **A result of the package belongs to an installation of R, not to the package alone.** The
  R of 2014 and the R of today draw the bootstrap samples differently, round differently, and
  on Intel processors add up in a higher precision than on Apple Silicon. On the package's own
  data the two installations agree when they are given the same samples. On small matrices
  full of ties they do not: of 20,000 random cases, 131 end with other genes in the network.
  So the port takes the installation as a parameter (`search.r2014.platform`); two are
  available, the two it was compared with.
- **The regression functions (`nodes.glm`, `nodes.fastLm`) are reproduced in their logic and
  to rounding error in their numbers**, not bit for bit: 2,192 of 2,198 cases agree to nine
  significant digits, the other six to six digits.
- **Four things cannot be reproduced, by anyone.** The figures. A run without
  `sampling.seed`. The neighbourhood of a seed that has the same value in every sample under
  a regression function, which R itself changes when the samples are listed in another order.
  And a user-written association function, which is R code.
- **Thirteen defects of the package are now known**, five of them found by this study. All
  are reproduced by the port and listed [below](#defects-of-the-package).

## What "the same" means

Three levels, stated for every comparison.

1. **Bit for bit.** Two numbers are the same when they are the same double: all 64 bits, with
   R's `NA` kept apart from `NaN`. Names are the same when they are the same strings in the
   same order, with a missing name where R has one. This is the level of the whole correlation
   path. No tolerance is involved.
2. **Character for character.** The files the package writes are compared after
   decompression, as text. (The compressed bytes carry a timestamp and differ between any two
   runs, of R as well.)
3. **To a stated tolerance.** Coefficients and p-values of the regression functions come out
   of different linear-algebra and distribution routines in the two languages. They are
   compared relative to their size with a tolerance of 1e-9. Everything discrete in those runs
   (which genes are in which neighbourhood, D, O, the seeds added and removed) is still
   compared exactly.

When R stops with an error, the port must stop too, with the text of R's error message.

## The verbatim port

`search.r2014` is a transliteration, function by function. The R statement each block stands
for is quoted in a comment next to it, and the places where the package misbehaves are marked
`R defect` with the number they have on this page.

| R | `search.r2014` |
|---|---|
| `SEARCH` | `core.SEARCH` |
| `SEARCH.randomise`, `.SEARCH.randomise.inner` | `randomise.SEARCH_randomise`, `SEARCH_randomise_inner` |
| `summarise.iterations` | `summarise.summarise_iterations` |
| `nodes.correlation` | `association.nodes_correlation` |
| `nodes.glm`, `nodes.fastLm` | `glm.nodes_glm`, `glm.nodes_fastLm` |
| `getO`, `getOmax`, `get.matrix.format`, `estimate.nodes.degree` | `core.getO`, `getOmax`, `get_matrix_format`, `estimate_nodes_degree` |
| `write.data.files` | `randomise.write_data_files` |
| `removeNode`, `removeByName` | `rsem.removeNode`, `remove_by_name_vector`, `remove_by_name_list` |
| `make.default.plots`, `make.networkplot`, `make.scatter.plot`, `get.adjacency.matrices`, `get.program.defaults` | not ported |

Arguments keep their R names with the dots written as underscores. An example, the first one
of the package manual:

```python
from search.r2014 import SEARCH_randomise, read_matrix, summarise_iterations

x = read_matrix("tests/data/breast/Breastdata1/mRNA_abundance.txt.gz")
SEARCH_randomise(
    data=x,
    seed_nodes=x.rownames[:10],
    universe_nodes=x.rownames,
    D_lower=0,
    D_higher=0,
    O_lower=0.1,
    max_iter=3,
    sd_threshold=0,
    n_boot=2,
    sampling_seed=1,
    output_directory="out",
    correlation_threshold=0.5,
    correlation_method="spearman",
)
```

It writes `out/output/connectivity__iteration_1_1.txt.gz` and the other files of the package,
with the content R writes.

Most defects of the package are not arithmetic. They come from what R does with a named
vector when a name is missing, with a condition that is `NA`, with a list element set to
`NULL`, with a table of one row. So the port rests on a few modules that reproduce R itself:

| module | what of R it reproduces |
|---|---|
| `rsem` | named vectors and lists, missing names, indexing by name and by condition |
| `rstats` | `sd`, `cor` (Pearson, Spearman, Kendall, every `use`), `round`, `p.adjust` |
| `rng` | `set.seed` and `sample`, as before and as after R 3.6.0 |
| `rformat` | how `write.table` writes a number and how `read.table` reads one |
| `rio` | `write.table`, the `gzip` call, `read.table`, as the package calls them |
| `x87` | the 80-bit arithmetic of Intel processors |
| `platform` | which installation of R all of the above stand for |

## How it was checked

**The package runs unmodified.** `validation/r_parity/r/r2014_oracle.R` reads the sources
under `reference/r_sources/SEARCH/R` as they are, runs a file of cases and writes everything
the package returns or the error it stops with. Two things are replaced: `mclapply` of the
removed `multicore` package by the one of `parallel`, and the plotting function by an empty
one. The script needs no R package, so it also runs in R 3.1.0.

**Cases.** Two kinds, 25,050 in all.

- *Named cases* run the package on Breastdata1, the data it ships with (1231 probesets, 77
  samples): the two examples of the manual as written there, one case for every defect,
  bootstrap and jackknife runs followed by `summarise.iterations`, the three correlation
  methods, the six GLM families, covariates and interactions, and calls that R refuses. Two
  more run on the GBM IT cohort and six on a small matrix with gene names that look like
  numbers.
- *Random cases* are small matrices (4 to 14 genes, 3 to 18 samples) with random arguments,
  built to reach the states in which the package misbehaves. With a handful of samples many
  correlations sit exactly on the threshold; missing values and constant rows are frequent;
  the universe is given in an order that is not the one of the matrix. One case in four ends
  with R stopping on an error.

**The building blocks, separately.** `r/r2014_primitives.R` asks R directly for random draws,
numbers read and written, rounding, standard deviations, correlations, adjusted p-values and
row means; `py/check_r2014_primitives.py` compares each with the port's version.

**Every function, separately.** `r/r2014_units.R` calls one function of the package per case,
with arguments chosen for that function: the ordinary ones and those at the edges (nothing,
one element, no names, missing names, a name given twice, missing values, ties).
`py/r2014_units.py` makes the cases, calls the port's function and compares in the same exact
sense. The rules of R the port wrote down (what `v[names] <- values` does with a name that is
not there, how `write.table` writes a number) are compared the same way, each as the R
expression it stands for. [Tests of the verbatim port](r2014_tests.md) lists the cases.

**The same in two installations of R.** Every script was run with R 4.4.3 on macOS arm64. The
building blocks and the named cases on Breastdata1 with correlations were also run with R
3.1.0 on Linux x86-64, in a container.

## Results

### R 4.4.3 on macOS arm64

"Equal" counts the cases in which nothing differs; "R stops" those in which R ends with an
error and the port ends with the same error.

| family | cases | equal | of which R stops | differ | numbers compared | largest difference |
|---|---|---|---|---|---|---|
| `SEARCH()` on Breastdata1 | 19 | 19 | 6 | 0 | 63,204 | 0 |
| `SEARCH()` on the GBM IT cohort | 2 | 2 | 0 | 0 | 35,984,987 | 0 |
| `SEARCH()` on random matrices | 20,000 | 20,000 | 4,996 | 0 | 1,047,655 | 0 |
| `SEARCH.randomise()` and the summary on Breastdata1 | 8 | 8 | 0 | 0 | 371,002 | 0 |
| `SEARCH.randomise()` and the summary on random matrices | 2,500 | 2,500 | 905 | 0 | 427,380 | 0 |
| `SEARCH.randomise()`, gene names that look like numbers | 6 | 6 | 0 | 0 | 1,548 | 0 |
| regression functions on Breastdata1 | 15 | 15 | 2 | 0 | 2,439,110 | 9.7e-12 |
| regression functions on random matrices | 2,183 | 2,177 | 185 | 6 | 633,602 | 7.9e-10 |
| regression functions, a constant seed | 317 | 88 | 40 | 229 | | see [R11](#defects-of-the-package) |

For the files, "numbers compared" counts lines. On the GBM IT cohort the port takes 9 seconds
for the three growth iterations and the package four minutes; both add 912, 819 and 297 seeds
and remove 106.

The six regression cases that differ are all `nodes.fastLm` with a model that has an
interaction between two covariates. Memberships, D and O are equal; what differs is the sixth
to ninth significant digit of some standard errors. `fastLm` takes them from the
pseudo-inverse of X'X, which squares the conditioning of the model, and the two languages
compute that pseudo-inverse with different routines.

The building blocks, each compared bit for bit or character for character:

| block | compared | different |
|---|---|---|
| `set.seed` and `sample`, both kinds: 9 seeds by 10 sizes each | 180 | 0 |
| a number as `write.table` writes it | 78,334 | 0 |
| a number as `read.table` reads it, 3 to 20 significant digits | 21,000 | 0 |
| `round(x, 5)` and `round(x, 3)` | 256,722 | 0 |
| `cor`: Pearson, Spearman, Kendall with each of the five `use` | 12,726 | 0 |
| `sd` | 7,590 | 0 |
| `p.adjust`, eight methods | 9,280 | 0 |
| `rowMeans` of ten values | 4,000 | 0 |

### Function by function

Run on 2026-11-20 in R 4.4.3 on macOS arm64: one call of one function per case. The first set
is the one kept in the test suite; the second is the same generator with twenty times the
random cases.

| | cases | equal | of which R stops | refused by the port, as stated | differ |
|---|---|---|---|---|---|
| the 14 functions of the package, test set | 1,062 | 1,058 | 153 | 4 | 0 |
| the 18 rules of R, test set | 574 | 573 | 23 | 1 | 0 |
| the 14 functions of the package, larger set | 16,526 | 16,520 | 2,021 | 4 | 2 |
| the 18 rules of R, larger set | 9,124 | 9,123 | 179 | 1 | 0 |

The two cases that differ are `nodes.fastLm` with covariates, in standard errors that agree
with R to eight significant digits and not to nine: the difference of the six whole-run cases
above. The five refusals are inputs the port does not model, tested as refusals.
[Tests of the verbatim port](r2014_tests.md) has the counts per function and the list of the
cases.

These numbers are those of the port after seven corrections. The first run of the unit cases
found 37 differences in 1,631 cases where the whole runs had found none, and the larger set one
more: `write.table` writes `NaN` as `NA` and the port wrote `NaN`; `cor` loses the names of
the genes when exactly one sample is complete under a rank method; three error messages were
not R's; three rules of R the package never reaches and the reading of numbers below 1e-307
were wrong. The list, with
what R does in each, is in [Tests of the verbatim port](r2014_tests.md#what-the-unit-tests-found).
After the corrections the comparison of the whole runs was repeated against the stored
answers of R, for both installations: every report is identical to the byte, so the tables
above and below stand.

### R 3.1.0 on Linux x86-64

With the port set to `platform.R_3_1_X86`:

| family | cases | equal | of which R stops | differ | numbers compared |
|---|---|---|---|---|---|
| `SEARCH()` on Breastdata1 | 19 | 19 | 6 | 0 | 63,204 |
| `SEARCH.randomise()` and the summary on Breastdata1 | 7 | 7 | 0 | 0 | 265,325 |
| the building blocks above, in R 3.1.0 | 46 blocks | 46 | | 0 | 389,671 |

Not run in R 3.1.0: the random cases, the regression cases and the GBM cohort. R 3.1.0 exists
only for x86 and runs emulated on this machine, about ten times slower; the run of the random
cases was under way when the container engine restarted, and it was not started again. Until
it is, the 80-bit arithmetic of the port is verified on every building block and on real
data, and not on the small matrices full of ties where it matters most.

## A result belongs to an installation of R

The same package on the same data gives different numbers in two installations of R. This is
a property of R, measured here because it decides what "reproducing the 2014 results" can
mean.

What changes between R 3.1.0 on x86-64 and R 4.4.3 on arm64:

| what | R 3.1.0, x86-64 | R 4.4.3, arm64 | how often they differ |
|---|---|---|---|
| `sample()` | index = floor(n * u) | by rejection, since R 3.6.0 | every bootstrap and jackknife sample |
| sums inside `cor`, `sd`, `rowMeans`, `sum` | 80-bit, rounded once at the end | 64-bit, rounded at every addition | `rowMeans` of ten values: 1,722 of 4,000; `sd`: 1,515 of 7,560 |
| `cor`, Pearson | | | 396 of 1,010 values, in the last digit |
| `cor`, Spearman | | | 29 of 1,010 values, in the last digit |
| `round(x, 5)` | rint(x * 10^5) / 10^5 in 80 bits | the nearer of the two neighbouring decimals, since R 4.0.0 | 256 of 128,361 values |
| a number read from text | exact up to 19 digits, then rounded twice | exact up to 15 digits | 15 digits or fewer: 2 of 9,337; 17 digits: 866 of 2,952 |
| a number written with 15 digits | scaled in 80 bits | scaled in 64 bits | 171 of 78,334 |
| `cor` below -1 by rounding | left as it is | set to -1 | rare |

Two of these reach the results of the package. D is written to the files rounded to five
decimals, and 0.003125 is 0.00312 for one R and 0.00313 for the other. And
`summarise.iterations` compares the mean of D over the replicates with the null values: the
mean of ten copies of 0.1 is 0.1 in 80 bits and 0.09999999999999999 in 64, so a null value of
0.1 reaches the mean on one machine and not on the other, and the p-value changes.

How much this matters was measured in two ways.

*In R, on the package's data.* The 26 cases that both installations ran:

| | cases | identical | differ in numbers only | differ in results |
|---|---|---|---|---|
| `SEARCH()` on Breastdata1, no random draw | 19 | 18 | 1 (Pearson, last digits) | 0 |
| `SEARCH.randomise()` with the samples of R 3.1.0 drawn in both | 1 | 1 | 0 | 0 |
| `SEARCH.randomise()` with each R's own `sample()` | 6 | 0 | | 6: other samples |

*With the port, on the random cases.* The port run as each of the two installations, with the
same random draws. This is a prediction of the port, verified in R only as far as the tables
above go:

| | cases | identical | numbers only | other genes | one stops, the other does not |
|---|---|---|---|---|---|
| `SEARCH()` on random matrices | 20,000 | 17,887 | 1,972 | 131 | 10 |
| `SEARCH.randomise()` and the summary on random matrices | 2,500 | 2,331 | 162 (the summary: mean, p-value) | 7 (the files) | |

The random matrices are a worst case by design. With 77 samples and continuous values, as in
Breastdata1, a correlation does not sit on the threshold and the two installations agree.

What follows for a result obtained with the package in 2014: its bootstrap samples are those
of `sample()` before R 3.6.0, which an R of today gives only after
`RNGkind(sample.kind = "Rounding")`; and if it was computed on an Intel machine, its last
digits, and through them some p-values of the summary, are those of 80-bit sums, which an
Apple Silicon machine does not produce in R either. The port reproduces both through
`platform.R_3_1_X86`.

```python
from search.r2014 import platform

with platform.use(platform.R_3_1_X86):
    SEARCH_randomise(...)
    summarise_iterations(...)
```

A third installation matters and was not compared: R 4 on Linux x86-64, as on a cluster. It
has the `sample()` and the `round()` of today with the 80-bit sums of Intel. The port has no
platform for it and refuses to guess.

## Defects of the package

Every one is reproduced by the port, and none is in the `search2014` preset of the engine. R1
to R8 are described with their measured effect in [Porting from R](porting.md); R9 to R13 were
found by this study.

| | where | what happens |
|---|---|---|
| R1 | `SEARCH.R:26-27` | The standard deviations of the sd prefilter are computed in the row order of the matrix and named in the order of the universe, seeds first: the filter removes the wrong genes. |
| R2 | `getOmax.R:25` | A seed is compared only with the seeds after it in the list; with `keep.seeds = FALSE` the last one always gets Omax 0. |
| R3 | `SEARCH.R:124` | With `keep.seeds = FALSE` an old seed removed on Omax keeps its neighbourhood: D gets missing names and the next iteration stops. |
| R4 | `SEARCH.R:204` | A gene removed from the universe stays a member of the other neighbourhoods and keeps counting in every later O. |
| R5 | `nodes.glm.R:107`, `nodes.fastLm.R:93` | A neighbourhood of one member loses its name; the seed is not counted in its own D and gets D = -1 / (K - 1). |
| R6 | `nodes.glm.R:109` | A seed with no significant association vanishes from the result; when every seed does, `SEARCH` stops. |
| R7 | `nodes.correlation.R:31` | A missing correlation puts an entry named `NA` in the neighbourhood, and the run stops when the matrices are built. |
| R8 | `summarise.iterations.R:124` | A gene present in a null run and in no real run stops the function. |
| R9 | `summarise.iterations.R:42-55` | The columns of the summary are numbered 1 to the number of files found, and each is filled from the replicate with that number. When a replicate in the middle has no file for an iteration, because it stopped earlier, the replicates after it are read and never used. With replicates 1 and 3 present, the mean is replicate 1 alone. The same for the null runs. |
| R10 | `write.data.files.R` | The path of every file is pasted into a shell command, `gzip --force --best <path>`, without quotes. With a blank in the output directory `gzip` fails, R does not notice, the files stay uncompressed, and `summarise.iterations`, which looks for `.gz`, finds none and stops. |
| R11 | `nodes.glm.R`, `nodes.fastLm.R` | A seed with the same value in every sample is fitted exactly by the intercept. The residuals are rounding noise, and so are the p-values that decide its neighbourhood. For a seed constant at 8.7 and 200 unrelated genes, R 4.4.3 finds 109 of them significant at p < 0.01; with the samples listed in reverse order, 100, of which 56 are the same. For a seed constant at 8, none. R 3.1.0 on x86 finds none in either case. |
| R12 | `nodes.fastLm.R` | A gene with the same value in every sample makes the model rank-deficient; `fastLm` returns a coefficient and a p-value for it all the same, so it can enter a neighbourhood and become a seed, and R11 then applies to it. (`nodes.glm` gives it no coefficient.) |
| R13 | `summarise.iterations.R:27` | `read.table` converts the column of gene names before using it as row names. When every name of a file looks like a number, the names come back as numbers: `007` as `7`, `1.50` as `1.5`, `100000` as `1e+05` if another name has a decimal point. Two names that become equal, or a gene called `NA`, stop the function. Names that are all `T`, `F`, `TRUE` or `FALSE` come back as `TRUE` and `FALSE`. |

Smaller points, also reproduced:

- `summarise.iterations` stops when no replicate reached the iteration asked for, when
  there is exactly one null file, when there is one gene, and with its own default arguments
  (`batches.NULL = NULL` makes `1:batches.NULL` an error). It cannot be called without null
  runs.
- `SEARCH.randomise(n.boot = 0)` runs two replicates, numbered 1 and 0: `1:0` counts down.
- With one seed, D of the seed is 0 / 0.
- A correlation matrix with a single column is written as a single row.
- O and the correlation matrix exist only for iterations followed by another one.
- The table of powers of ten that R 3.1.0 uses to print numbers on x86 is declared exact up
  to 10^27 and is not from 10^23 on. This is a defect of R, not of the package, and it had to
  be reproduced to print as that R prints.

## What is not reproduced

- **The figures.** The plotting functions need the igraph of 2014 and are not ported. In the
  package they run after the files of a replicate are written; if one of them failed, the run
  stopped there, and the port does not know.
- **Runs without `sampling.seed`.** R then continues the random stream of the session.
  Nothing can reproduce that, R included.
- **A constant seed with a regression function** (R11). The port computes a neighbourhood of
  the same kind and raises the warning `NotReproducible`.
- **The regression functions bit for bit.** They agree with R to nine digits, six for the
  standard errors of badly conditioned models. A p-value that close to the threshold could
  fall on the other side; it did not in 2,198 cases.
- **A seed in its own neighbourhood, under a regression function with a threshold below about
  1e-140.** A seed is regressed on itself like on any other gene; the fit is exact and its
  p-value is rounding noise: 0 in R for the two seeds of a test case, between 8e-147 and
  4e-173 in the port. With the thresholds in use (0.01, 0.05) the seed always passes. With
  1e-200, R kept both seeds and the port none or one.
- **`formula.string` in general.** It is R code. The port reads sums of covariate names and
  their interactions (`a + b`, `a:b`, `a*b`, also with `x.node`) and refuses the rest.
  Covariates that are text become factors with their levels in the order of the character
  codes, which is R's order in the C locale and not in every other.
- **Association functions written by the user.** `nodes.association.function` can name any R
  function; the port has the three of the package and accepts a Python function.
- **`mc.cores` above 1 when a fit fails.** With one core R stops; with more it goes on with
  whatever the failed worker returned. The port behaves as with one core.
- **States that need inputs the package was not meant for.** Duplicated row names or seeds, a
  file with quotes or `#` in a gene name, an output path with characters the shell
  interprets other than blanks. The port stops with `UnportedState` instead of guessing. It
  did not happen in any of the 25,050 cases.
- **What R prints on the console.**

## Not compared

- R 3.1.0 on the random cases, the regression functions and the GBM cohort (see above).
- R 3.1.0 function by function: the unit cases were answered by R 4.4.3 only.
- R 4 on Linux x86-64 and any R on Windows.
- Versions of R between 3.1.0 and 4.4.3. Two changes are dated in R's own news (`sample()` in
  3.6.0, `round()` in 4.0.0); when `cor` began to clamp at -1 was not looked up.
- Data files with more than 15 significant digits. The reader of the port reproduces the one
  of R for them too (21,000 numbers of 3 to 20 digits), but the matrices of this project have at most 15.

## Repeating the study

```bash
export R_LIBS=/path/to/a/library/with/RcppArmadillo      # optional: the nodes.fastLm cases
export SEARCH_GBM_DATA=/path/to/GBM_HYPOXIA_SIGNATURE     # optional: the two GBM cases
validation/r_parity/run_r2014.sh            # 25,050 cases in the R of the machine, 15 minutes
validation/r_parity/run_r2014_r31.sh 2      # the same cases in R 3.1.0, in containers
validation/r_parity/run_r2014_units.sh 20   # every function on its own, R of the machine

Rscript validation/r_parity/r/r2014_primitives.R out/primitives
.venv/bin/python validation/r_parity/py/check_r2014_primitives.py out/primitives/primitives.json
.venv/bin/python validation/r_parity/py/compare_r_versions.py validation/r_parity/out/r2014
.venv/bin/python validation/r_parity/py/compare_platforms.py validation/r_parity/out/r2014
```

The reports are under `validation/r_parity/results`: `report_r2014.md` and
`report_r2014_r31.md` with one row per family, `report_r2014_cases.csv` with the named cases
and every random case that is not simply equal, `report_r2014_primitives*.csv` for the
building blocks. The test suite does not need R: `tests/unit/test_r2014_golden.py` compares
the port with answers recorded from both installations (436 cases from R 4.4.3, 9 from R
3.1.0, and a sample of every building block), and `tests/unit/test_r2014_defects.py` works
each defect out on an example small enough to follow by hand.
