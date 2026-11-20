# Verbatim port against R version 3.1.0 (2014-04-10) (x86_64-unknown-linux-gnu)

Port run as: R 3.1.0, Linux x86-64. R side: `sample()` kind Rounding, long double of 16 bytes, collation en_US.UTF-8, run on 2026-11-19.

| family | cases | equal | equal, R stops | differs | not ported | quantities | values | largest difference among the equal |
|---|---|---|---|---|---|---|---|---|
| SEARCH(), Breastdata1 | 19 | 13 | 6 | 0 | 0 | 510 | 63204 | 0 |
| SEARCH.randomise(), Breastdata1 | 7 | 7 | 0 | 0 | 0 | 570 | 265325 | 0 |
