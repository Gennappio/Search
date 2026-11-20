# How the verbatim port is tested

An introduction for readers who know what the SEARCH package does and have not followed the
port. The reference, with every case, table and finding, is
[Tests of the verbatim port](r2014_tests.md).

## What is being tested

The SEARCH R package (Haider and Buffa, version 0.0.2, 2014) was rewritten in Python twice.

- The **engine** of Search computes the method as it was meant, with the defects of the
  package corrected.
- The **verbatim port**, `search.r2014`, reproduces the package as it was written, defects
  included. Its purpose is to obtain again, without R, exactly what the 2014 package gives,
  and to tell whether a result depends on one of its defects.

These tests concern the verbatim port only. They answer one question: **does each function of
the port return exactly what the same function of the 2014 package returns in R, given the
same input?** For this code R is the specification. A test passes when the port gives R's
answer, also where R's answer is wrong.

## One test, from start to finish

`getOmax` receives a list of neighbourhoods (for each seed gene, the genes associated with
it) and gives each its largest overlap with the others: 1 when they have the same genes, 0
when they have none in common.

Given three identical neighbourhoods, the package in R answers:

```r
x <- list(s1 = c(g1 = 1, g2 = 1), s2 = c(g1 = 1, g2 = 1), s3 = c(g1 = 1, g2 = 1))
getOmax(x, x)
# 1 1 0
```

The right answer would be `1 1 1`. The package compares each neighbourhood only with those
after it, so the last one is compared with nothing and gets 0 (defect R2 of the package). The
port answers:

```python
getOmax(x, x).values
# [1.0, 1.0, 0.0]
```

It is wrong in the same way, so the test passes. If someone corrected the defect in the port,
this test would fail.

## How the tests work

1. **Cases.** A case is one function of the package and one input for it. There are 1,636.
2. **R answers once.** The unmodified R code of 2014 was run on every case, and its answers
   were saved in the repository, in `tests/data/r2014/units_r44_arm64.json.gz`.
3. **Python is compared every time.** The test suite calls the port's function on the same
   input and checks that the answer is identical to the one R gave.

R is needed only for step 2. Running the tests needs Python alone, so they run anywhere,
including the repository's continuous integration, in about 15 seconds. To confirm that the
saved answers really are R's, one script runs R again on every case and compares the port
with the fresh answers (see [Running them](#running-them)).

Numbers travel to R as hexadecimal floating point (`0x1.0000000000000p+0` is 1), which both
languages read exactly, and R's answers are saved with 17 significant digits, which is enough
to read back the same number. A test checks that a number sent to R comes back unchanged.

## What "identical" means

- **Numbers**: the same double, all 64 bits. R's missing value `NA` and `NaN` are different;
  `0` and `-0` are different.
- **Names**: the same names in the same order; missing where R has them missing; absent where
  R has none.
- **Structure**: the same kind of object (nothing, vector, matrix, list) with the same
  dimensions.
- **Files**: when a function writes files, the same files with the same text.
- **Errors**: where R stops with an error, the port stops with the same message. This happens
  in 176 of the cases.

Two exceptions. The coefficients and p-values of the two regression functions, `nodes.glm`
and `nodes.fastLm`, are compared to six significant digits, because the two languages compute
them with different linear-algebra routines; which genes they select is still compared
exactly. And the type of an empty vector is not compared: R tells an empty logical vector
from an empty character vector, Python does not.

## What is tested

### The 14 functions of the package (1,062 cases)

| R function | what it does | cases |
|---|---|---|
| `SEARCH` | the whole method: starts from the seeds, finds the associated genes, adds and removes seeds, repeats | 78 |
| `SEARCH.randomise` | runs `SEARCH` on resampled data (bootstrap) and writes the results to files | 36 |
| `.SEARCH.randomise.inner` | one replicate of the above | 22 |
| `summarise.iterations` | reads those files back and computes means and p-values | 103 |
| `write.data.files` | writes the result files | 81 |
| `nodes.correlation` | for each seed, the genes whose correlation passes the threshold: its neighbourhood | 189 |
| `nodes.glm`, `nodes.fastLm` | the same with a regression instead of a correlation | 70, 48 |
| `getO` | the overlap of two neighbourhoods | 73 |
| `getOmax` | the example above | 90 |
| `get.matrix.format` | builds the overlap and correlation matrices | 99 |
| `estimate.nodes.degree` | for each seed, how many other seeds its neighbourhood holds | 60 |
| `removeNode`, `removeByName` | remove elements from a list | 49, 64 |

Each function gets ordinary inputs and awkward ones: empty, a single element, no names,
missing or repeated names, missing values, ties. That is where R and Python part ways.

The package has 19 functions; the other five draw figures and were not ported. A test reads
the R sources and fails if a function is neither tested nor listed as not ported.

### The rules of R behind the defects (574 cases)

Most defects of the package come from what R does in unusual situations: a name that is not
there, a condition that is missing, an empty selection. Python behaves differently in those
situations, so the port had to write R's behaviour down, and each rule is tested on its own.
For example:

```r
v <- c(g1 = 1, g2 = 2, g3 = 3)
v[c(NA, NA, "q", "q", "g1")] <- c(10, 11, 12, 13, 14)
# g1 = 14, g2 = 2, g3 = 3, NA = 10, NA = 11, q = 13
```

R appends one element for each missing name and one for the new name `q` given twice
(keeping its last value), and updates `g1`. Eighteen rules are tested like this, among them
how R writes numbers to a file and reads them back.

### Tests about the tests

- **Every defect is covered.** For 12 of the 13 known defects of the package there are cases
  written for it, and a test checks that the defect is visible in R's answer and that the port
  reproduces it. The thirteenth (R11) has no case, because R's own answer there is rounding
  noise.
- **A repair is noticed.** Seven tests deliberately correct a defect inside the port, or
  replace a rule of R with the reasonable behaviour, and require that some case fails. This
  shows that the cases hold the port to what R does, not to what would be right.
- **The cases are complete and untouched.** Every function has enough cases, including cases
  in which it stops with an error, and the stored cases are exactly those the generator
  produces.
- **The comparison works.** It must report a difference of one bit in one number, one changed
  name, `NaN` in place of `NA`, and an error in place of a value.

## Where the cases come from

The package has no tests of its own, and the manual pages of its single functions have empty
examples. So the inputs were written for this purpose. Only the inputs: the expected answer of
every case was given by R.

| kind | cases | how they were made |
|---|---|---|
| chosen | 364 | written one by one after reading the R code of each function, each with a note on what it tries; at least one for every known defect |
| searched | 6 | for two defects, random inputs were tried until the defect appeared |
| random | 1,266 | generated by a program from a fixed seed, so always the same: small matrices (4 to 14 genes, 3 to 18 samples) with missing values, constant rows and many ties. No patient data |

The number of random cases per function was a practical choice: enough to reach the edges,
few enough for the tests to run in seconds. To check that they are enough, the same generator
was run with twenty times the random cases (25,650 cases, two hours of R). That run found one
difference the smaller set had not, an error message worded differently. It was corrected and
a case was added for it.

## Results

Last run on 2026-11-20, Python 3.12.9 and NumPy 2.5.3 on macOS (arm64), commit `276f700`:

| what | result |
|---|---|
| the tests of the verbatim port | 1,688 passed, none failed |
| the whole test suite of the project | 1,784 passed; continuous integration green on Linux and macOS, Python 3.12 and 3.13 |
| the 1,636 cases, R 4.4.3 run again | all as expected: 1,455 equal, 176 equal where R stops, 5 inputs the port refuses, as the tests require |
| 25,650 cases, R 4.4.3 | all as expected except 2 `nodes.fastLm` cases, whose standard errors agree to 8 significant digits instead of 9 (6 are required) |

These tests found seven places where the port did not do what R does, which 22,535 complete
runs of `SEARCH()` and `SEARCH.randomise()` had not reached. The most concrete: R writes a
`NaN` to a file as `NA`, and the port wrote `NaN`. All seven were corrected, and the complete
runs were compared with R's saved answers again afterwards, with no change.

## What the tests do not cover

- **The R of 2014.** The answers come from R 4.4.3 on macOS (arm64). R 3.1.0, released with
  the package, draws random samples and rounds differently, and on Intel processors adds in a
  higher precision. The port can imitate it, but function by function it was compared with R
  4.4.3 only; with R 3.1.0 it was compared on complete runs on the package's own data and on
  the basic operations.
- **Regression functions to the last bit**, as said above.
- **A seed regressed on itself with a significance threshold below about 1e-140.** Its p-value
  is rounding noise, and there the port and R disagree. With the usual thresholds (0.01,
  0.05) they agree.
- **The figures** of the package, which were not ported, and what R prints on the console.
- **Five unusual inputs**, such as two genes with the same name, which the port refuses
  instead of guessing. The tests check that it refuses them.

## Running them

```bash
# the tests: Python only, about 15 seconds
.venv/bin/python -m pytest tests/unit/test_r2014_units.py tests/unit/test_r2014_golden.py \
    tests/unit/test_r2014_defects.py

# the tests of one function
.venv/bin/python -m pytest tests/unit/test_r2014_units.py -k getOmax -v

# confirm the saved answers against R, one minute: needs R with the package jsonlite,
# and RcppArmadillo for the nodes.fastLm cases
validation/r_parity/run_r2014_units.sh
```

## Where things are

| file | what it is |
|---|---|
| `src/search/r2014/` | the verbatim port |
| `tests/unit/test_r2014_units.py` | the tests described here |
| `tests/data/r2014/units_r44_arm64.json.gz` | the 1,636 cases with R's answers |
| `tests/unit/test_r2014_golden.py` | complete runs and basic operations against R 4.4.3 and R 3.1.0 |
| `tests/unit/test_r2014_defects.py` | the defects and rules of R on examples worked out by hand |
| `validation/r_parity/py/r2014_units.py` | the generator of the cases and the comparison |
| `validation/r_parity/r/r2014_units.R` | the R side: one function of the package per case |
| `validation/r_parity/results/report_r2014_units.md`, `.csv` | the last comparison with R, per function and per case |
| [Tests of the verbatim port](r2014_tests.md) | the reference: every table, every defect case, every finding |
| [Replicating the 2014 package](replication.md) | the study of whether a result of the 2014 package can be reproduced |
