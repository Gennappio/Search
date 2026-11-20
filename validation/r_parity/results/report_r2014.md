# Verbatim port against R version 4.4.3 (2025-02-28) (aarch64-apple-darwin20)

Port run as: R 4.4.3, macOS arm64. R side: `sample()` kind Rejection, long double of 8 bytes, collation C.UTF-8, run on 2026-11-19.

| family | cases | equal | equal, R stops | differs | not ported | quantities | values | largest difference among the equal |
|---|---|---|---|---|---|---|---|---|
| SEARCH(), Breastdata1 | 19 | 13 | 6 | 0 | 0 | 510 | 63204 | 0 |
| SEARCH(), GBM IT cohort | 2 | 2 | 0 | 0 | 0 | 3358 | 35984987 | 0 |
| SEARCH(), random matrices | 20000 | 15004 | 4996 | 0 | 0 | 329722 | 1047655 | 0 |
| SEARCH.randomise(), Breastdata1 | 8 | 8 | 0 | 0 | 0 | 704 | 371002 | 0 |
| SEARCH.randomise(), gene names that look like numbers | 6 | 6 | 0 | 0 | 0 | 270 | 1548 | 0 |
| SEARCH.randomise(), random matrices | 2500 | 1595 | 905 | 0 | 0 | 74524 | 427380 | 0 |
| regression, Breastdata1 | 15 | 13 | 2 | 0 | 0 | 607 | 2439110 | 9.66e-12 |
| regression, constant seed | 317 | 48 | 40 | 229 | 0 | 7113 | 109070 | 1.74e-10 |
| regression, random matrices | 2183 | 1992 | 185 | 6 | 0 | 52299 | 633602 | 7.94e-10 |
