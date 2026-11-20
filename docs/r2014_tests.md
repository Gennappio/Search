# Tests of the verbatim port

`search.r2014` is the SEARCH 0.0.2 package of 2014 written again in Python line by line, with
its defects. For this code R is the specification: a test passes when the port returns what the
R package returns, not when it returns what would be right. This page says which tests exist,
what each one compares, how to run them and what the last run gave.
[Replicating the 2014 package](replication.md) has the study the tests come from.

## The three test files

| file | what it holds the port to | tests | needs R |
|---|---|---|---|
| `tests/unit/test_r2014_units.py` | every function of the package, and every rule of R the port had to write down, called on its own: 1,636 calls with the answer R gave to each | 1,664 | no |
| `tests/unit/test_r2014_golden.py` | whole runs of `SEARCH`, `SEARCH.randomise` and `summarise.iterations` (445 cases) and the functions of R the port reimplements (`sample`, reading and writing numbers, `round`, `sd`, `cor`, `p.adjust`, `rowMeans`), in R 4.4.3 on arm64 and in R 3.1.0 on x86-64 | 4 | no |
| `tests/unit/test_r2014_defects.py` | the rules of R and the defects of the package, each on an example small enough to work out by hand: these say what is expected and why | 20 | no |

None of them needs R: the answers were recorded from R once and are stored under
`tests/data/r2014`. R is needed only to record them again.

## Running them

```bash
.venv/bin/python -m pytest tests/unit/test_r2014_units.py tests/unit/test_r2014_golden.py \
    tests/unit/test_r2014_defects.py                   # all of them, 15 seconds
.venv/bin/python -m pytest tests/unit/test_r2014_units.py -k getOmax -v     # one function
.venv/bin/python -m pytest tests/unit/test_r2014_units.py -k "summarise_iterations/002" -v
```

Last run, 2026-11-20, Python 3.12.9 and NumPy 2.5.3 on macOS arm64: 1,688 passed, none
failed, none skipped. They are part of the whole suite (`pytest`: 1,784 passed without the GBM
data), which CI runs on Linux and macOS.

A failing unit test names the function, the case and the first thing that differs:

```
test_the_function_returns_what_it_returns_in_r[write_data_files/016]
write.data.files (an O of NaN): file shared_neighbourhood__iteration_1_1.txt.gz:
line 2: R 'g1\tNA' / port 'g1\tNaN'
```

## The unit tests

### What a case is

A case names one function and gives its arguments. `validation/r_parity/r/r2014_units.R`
reads the sources of the package as they are (`reference/r_sources/SEARCH/R`), calls that
function with those arguments and writes down what it returns, the files it writes, or the
message of the error it stops with. The test calls the port's function with the same
arguments and compares.

Two things of the package are replaced when R runs it, as everywhere in this project:
`mclapply` of the removed `multicore` package by the one of `parallel`, and the plotting
function by an empty one.

The numbers of a case reach R as hexadecimal floating point, which R reads exactly, and R's
answers are stored with 17 significant digits, which is enough to read back the same double.
A test checks this path: what R returns for the rule `v` is what was sent, bit for bit,
including 5e-324, the largest double, -0, `NA`, `NaN` and both infinities.

### What "the same" means

- **Names**: the same strings in the same order; missing (`NA`) where R has them missing;
  absent (`NULL`) where R has none. A vector without names and a vector whose names are all
  missing are different things, in R and in the test.
- **Numbers**: the same double, all 64 bits. `NA` and `NaN` are different, -0 and 0 are
  different.
- **Shape**: `NULL`, a vector, a matrix with its dimensions and both sets of names, a list
  with its names, compared recursively.
- **Files**: the same file names and, once uncompressed, the same characters.
- **Errors**: where R stops, the port stops with the same message.

Two exceptions, both stated in the code of the comparison (`r2014_units.differences`):

- The type of a vector without elements is not compared. R tells `logical(0)` from
  `character(0)` and from `list()`; the port keeps one empty Python list for all three.
  Whether the empty vector has names is compared.
- The coefficients and p-values of `nodes.glm` and `nodes.fastLm` are compared to six
  significant digits (nine when the cases were recorded, on one machine). They come out of
  other linear-algebra routines than R's. Their names, shapes and memberships are exact.
  Where a seed is regressed on itself the fit is exact and its standard error is rounding
  noise on both sides: there only the estimate is compared.

Three tests check the comparison itself: it reports a result that differs from R's in the last
bit of one number, in one name, in a missing name, in `NaN` for a number, in names removed,
and in an error for a value; it tells `NA` from `NaN` and -0 from 0.

### The functions of the package

The package defines 19 functions. Fourteen are ported and have cases; five draw the figures
and are not ported (`make.default.plots`, `make.networkplot`, `make.scatter.plot`,
`get.adjacency.matrices`, `get.program.defaults`). A test reads the R sources and fails if a
function is in neither list.

"Chosen" cases are written one by one for a situation, each with a note that says which.
"Random" cases are generated from a fixed random seed.

| R function | port | cases (chosen + random) | R stops in | what the chosen cases cover |
|---|---|---|---|---|
| `getO` | `core.getO` | 13 + 60 | 0 | no member, all members, a member twice, no names, missing names, empty on one or both sides (0 / 0) |
| `getOmax` | `core.getOmax` | 10 + 80 | 0 | no list (`NULL`), one list, identical lists, x longer than y, lists without names |
| `removeNode` | `rsem.removeNode` | 9 + 40 | 0 | nothing to remove, `NULL`, everything, an element twice, a missing value |
| `removeByName` | `rsem.remove_by_name_vector`, `remove_by_name_list` | 14 + 50 | 0 | a vector and a list; no names; a name twice; a missing name |
| `get.matrix.format` | `core.get_matrix_format` | 19 + 80 | 11 | one list; one member (names dropped); a column of missing values (dropped); a member that is not a feature; a member twice; empty lists |
| `estimate.nodes.degree` | `core.estimate_nodes_degree` | 10 + 50 | 0 | one seed among the members; none; a matrix without names; no column; `NaN` |
| `nodes.correlation` | `association.nodes_correlation` | 39 + 150 | 26 | the three methods and the five rules for missing values, every combination with a missing value; the universe in another order, smaller, larger; a seed outside it; one gene; arguments that do not exist |
| `nodes.glm` | `glm.nodes_glm` | 30 + 40 | 5 | the six families; the eight p-value adjustments; six formulas with covariates and interactions; a neighbourhood of one and of none; a model with one coefficient fewer for one gene |
| `nodes.fastLm` | `glm.nodes_fastLm` | 24 + 24 | 4 | the same without the families; a gene with one value in every sample |
| `SEARCH` | `core.SEARCH` | 12 + 66 | 29 | the sd prefilter; one seed; no seed; `max.iter = 0`; three iterations; a missing value; cases searched for that show R3 and R4 |
| `write.data.files` | `randomise.write_data_files` | 21 + 60 | 18 | rounding and sorting of D with ties, `NA`, `NaN`, `Inf`; matrices of one row, one column, no column; names missing or twice; coefficients of 15 digits |
| `.SEARCH.randomise.inner` | `randomise.SEARCH_randomise_inner` | 6 + 16 | 7 | bootstrap, jackknife, a sampling that does not exist, the seed as text, replicate 0 |
| `SEARCH.randomise` | `randomise.SEARCH_randomise` | 12 + 24 | 13 | a matrix of one row or column; no seed in the matrix; seeds twice; `n.boot = 0`; replicates 4 to 6; a blank in the output directory |
| `summarise.iterations` | `summarise.summarise_iterations` | 23 + 80 | 40 | a replicate without a file; a gene only in a null run; one null file; one gene; no null run; two iterations; equal means; ten values of 0.1; gene names that look like numbers |
| | | **1,062** | **153** | |

The random cases of `nodes.correlation` go through the 15 combinations of method and rule for
missing values in turn, on matrices with missing values in three cases out of ten. Those of
`summarise.iterations` are sets of connectivity files written by hand, with replicates
missing, genes that differ between replicates and values that tie with the null ones.

### The rules of R

Most defects of the package come from what R does at the edges: a name that is not there, a
condition that is missing, an empty selection. The port wrote each of those rules down once,
and each is tested as the R expression it stands for.

| R expression | port | cases | R stops in |
|---|---|---|---|
| `v[names] <- values` | `rsem.set_by_names` | 36 | 0 |
| `v[names]` | `rsem.get_by_names` | 34 | 0 |
| `v[v < t]`, `v[v > t]`, `v[abs(v) > t]` | `rsem.subset_logical`, `rsem.compare` | 111 | 0 |
| `names(v) <- names` | `rsem.set_names` | 38 | 16 |
| `union`, `intersect`, `setdiff`, `%in%`, `unique` | the same names in `rsem` | 172 | 0 |
| `c(a, b)` for two lists | `rsem.concat_lists` | 35 | 0 |
| `sort(v, decreasing = TRUE)` | `randomise._sort_decreasing` | 35 | 0 |
| `as.numeric(sprintf("%.3f", v))` | `randomise._three_decimals` | 5 (813 numbers) | 0 |
| `a:b` | `summarise._colon` | 7 | 1 |
| `write.table(x, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")` | `rio.table_of_vector`, `table_of_matrix`, `table_of_strings` | 51 | 1 |
| `as.matrix(read.table(file, row.names = 1, header = TRUE, sep = "\t"))` | `rio.read_connectivity` | 47 | 5 |
| `v` (the numbers travel exactly) | | 3 (315 numbers) | 0 |
| | | **574** | **23** |

`sample`, `round`, `sd`, `cor`, `p.adjust`, `rowMeans` and the reading and writing of single
numbers are not in this table: the study compared 389,832 values of them with R, and
`test_r2014_golden.py` keeps a sample of those for both installations of R.

### Every defect of the package, as R shows it

For each defect of the list in [Replicating the 2014 package](replication.md#defects-of-the-package)
there are cases written for it. A test checks two things on each: that the defect is visible
in what R returned (the check is in the last column), and that the port returns the same.

| defect | cases | what is checked in R's answer |
|---|---|---|
| R1 | `SEARCH/002` | the flat gene is in D, the seed with a large standard deviation is not |
| R2 | `getOmax/001`, `002` | identical lists: the first gets 1, the last gets 0 |
| R3 | `SEARCH/013` to `015` | a missing name among the seeds added |
| R4 | `SEARCH/016` to `018` | a removed gene is still a member of a neighbourhood |
| R5 | `nodes_glm/004`, `nodes_fastLm/004` | a neighbourhood of one value without a name |
| R6 | `nodes_glm/003`, `nodes_fastLm/003` | fewer neighbourhoods than seeds |
| R7 | `nodes_correlation/019` to `021`, `SEARCH/012` | a member with a missing name; `SEARCH` stops |
| R8 | `summarise_iterations/004` | the function stops |
| R9 | `summarise_iterations/002`, `003` | a column no replicate filled, with a file left unused |
| R10 | `SEARCH_randomise/011` | no file is compressed |
| R12 | `nodes_fastLm/023` | the constant gene has a coefficient |
| R13 | `summarise_iterations/016` to `020` | the genes of the summary are not called as in the files, or the function stops |

R11 has no case, and cannot have one: for a seed with the same value in every sample under a
regression function R's own answer is rounding noise, so there is nothing to record. The
port warns (`NotReproducible`), and `test_r2014_defects.py` tests that it does.

The cases for R3 and R4 were found by trying random inputs on the port until the symptom
appeared, and then confirmed in R.

### The tests hold the port to R, not to what is right

Seven tests repair something in the port and run the cases again. Each must make at least one
case fail; if none failed, the cases would not be testing that behaviour.

| what is repaired | cases that must notice |
|---|---|
| R2: `getOmax` compares every pair of lists | `getOmax` |
| R7: a missing correlation is left out of the neighbourhood | `nodes.correlation` |
| R9: every file found is a column of the summary | `summarise.iterations` |
| R10: the files are compressed whatever the path | `SEARCH.randomise` |
| R13: gene names are read as they are written | `summarise.iterations` |
| `write.table` writes `NaN` as `NaN` (R writes `NA`) | `write.data.files` |
| equal values are sorted in the opposite order | `write.data.files` |

### The stored cases are complete and are not edited by hand

- Every function of the package has cases or is in the list of those not ported (see above).
- Every rule the port's comparison knows has cases, and the other way round.
- Every function has at least a set number of cases, and every function that can stop has
  cases in which it does.
- The cases written one by one in the stored file are, one by one, those the generator
  (`validation/r_parity/py/r2014_units.py`) makes, and the file has as many random ones for
  every function.
- The answers come from the R the port stands for (R 4.4.3, `long double` of 8 bytes), and
  R skipped no case.

### Inputs the port refuses

For five of the 1,636 cases the test does not expect R's answer but a refusal: the port stops
with `UnportedState` instead of guessing. None of these states can be reached through
`SEARCH` or `SEARCH.randomise` with a data matrix whose gene names are all different.

| case | input | R |
|---|---|---|
| `get_matrix_format/016` | neighbourhood lists without names | stops |
| `get_matrix_format/017` | two lists with the same name | returns a matrix |
| `get_matrix_format/018` | a list with a missing name | stops |
| `get_matrix_format/019` | a gene twice among the row names | returns a matrix |
| `assign_by_name/006` | `v[names] <- values` on a vector without names | returns a vector |

## What the unit tests found

The whole-run comparison had found no difference in 22,535 cases. The first run of the unit
cases found 37 differences in 1,631. Fifteen came from one mistake in the new cases, two from
the finding described after the table, and twenty from six places where the port was not doing
what R does, none of which a whole run had reached. A larger run found a seventh (below). All
seven are corrected, and the whole-run comparison was repeated afterwards against the stored
answers of R, with reports identical to the byte.

| what | R | the port, before | when it matters |
|---|---|---|---|
| `write.table` and `NaN` | writes `NA`: its test for a missing value is `ISNAN` | wrote `NaN` | an O of two empty neighbourhoods (0 / 0), a `NaN` among the coefficients of a regression |
| `cor` with Spearman or Kendall, `complete.obs` or `na.or.complete`, and exactly one complete sample | returns a matrix without names (its `Rank` helper returns `row(u)` for one row), so the neighbourhood is missing values without names | missing values under missing names | data with missing values and few samples |
| `nodes.glm` and `nodes.fastLm` with a universe of one gene | stop at `rep(NA, ncol(x))`: "invalid 'times' argument" | stopped later, with "incorrect number of dimensions" | only the message |
| a correlation method that does not exist | the message quotes the choices with double quotes (`dQuote`) | single quotes | only the message |
| three rules outside what the package does: `v[names]` on a vector without names, `names(NULL) <- character(0)`, `c()` of two empty lists | no names; an error; no names | missing names; `NULL`; empty names | never through the package |
| reading a number below 1e-307 | scales it in two steps and reads 5e-324 as 5e-324 | read 0 | never in a file of the package |
| a regression in which one gene's model has fewer coefficients than the others (a covariate level seen only in a sample where that gene is missing) | fits every gene, then `cbind` stops and names the first table that does not fit: "(see arg 6)" | stopped at that gene and always said "(see arg 2)" | only the message |

One finding is not a difference of the port. With a regression function a seed is also
regressed on itself. The fit is exact, the residuals are rounding noise, and so is the
p-value: for the two seeds of one case R returned 0 and 0, the port 4e-173 and 8e-171 with
`nodes.glm`, 8e-147 and 0 with `nodes.fastLm`. With a `significance.threshold` of 1e-200 R
therefore kept both seeds in their own neighbourhoods and the port kept none, or one. With
any threshold in use (0.01, 0.05) the seed passes on both sides, and so it does at 1e-12, the
value of the case in the tests. A threshold below about 1e-140 is the limit: there the port
does not reproduce R.

### The larger run

The same generator with twenty times the random cases gives 25,650 cases. R answered them in
2 hours 14 minutes. That run found the last row of the table above, in one case of 829 of
`nodes.glm`. With it corrected:

| | cases | equal | of which R stops | refused, as stated | differ |
|---|---|---|---|---|---|
| the 14 functions of the package | 16,526 | 16,520 | 2,021 | 4 | 2 |
| the 18 rules of R | 9,124 | 9,123 | 179 | 1 | 0 |

The two cases that differ are `nodes.fastLm` with covariates: some standard errors agree with
R to eight significant digits and not to nine (largest relative difference 6.7e-9). This is
the difference already known from the whole runs, and it is inside the six digits the tests
ask of the regression functions: at that tolerance all 1,332 regression cases of the run are
equal. `validation/r_parity/results/report_r2014_units_x20.md` has the counts per function.

## What these tests do not cover

- **What R prints and warns.** Only what a function returns, writes or stops with is compared.
- **R 3.1.0 function by function.** The unit cases were recorded in R 4.4.3 on macOS arm64
  only. For R 3.1.0 on x86-64 the tests have the whole runs on Breastdata1 (9 cases) and the
  building blocks (389,671 values). Recording the unit cases there needs R 3.1.0 in a
  container.
- **The regression functions bit for bit**, as said above.
- **The figures**, which are not ported.
- **Runs without `sampling.seed`**, which R itself cannot repeat.
- **The GBM cohort.** Nothing derived from patient data is in the test data. The comparison on
  the GBM IT cohort is in the study and is run by hand.
- **Inputs the port refuses**, listed above and in
  [Replicating the 2014 package](replication.md#what-is-not-reproduced).

## Recording the answers again

```bash
export R_LIBS=/path/to/a/library/with/RcppArmadillo    # for the nodes.fastLm cases
validation/r_parity/run_r2014_units.sh        # the 1,636 cases, one minute; rewrites tests/data/r2014/units_r44_arm64.json.gz
validation/r_parity/run_r2014_units.sh 20     # twenty times the random cases, two hours; the test data are left alone

# the port changed and R's answers are still on disk: compare again, without R
.venv/bin/python validation/r_parity/py/r2014_units.py check validation/r_parity/out/r2014_units
```

It needs R with `jsonlite`. The script prints one line per function and every case that
differs, writes `validation/r_parity/results/report_r2014_units.md` with one row per function
and `report_r2014_units.csv` with one row per case, and ends with an error if any case
differs. The answers of the whole runs and of the building
blocks are recorded by `make_r2014_golden.py`, described in `validation/r_parity/README.md`.
