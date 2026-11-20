Unit cases of `search.r2014` against R version 4.4.3 (2025-02-28) on aarch64-apple-darwin20, 2026-11-20.

| function or rule of R | cases | equal | of which R stops | refused by the port, as stated | differ | skipped |
|---|---|---|---|---|---|---|
| `getO` | 1213 | 1213 | 0 | 0 | 0 | 0 |
| `getOmax` | 1610 | 1610 | 0 | 0 | 0 | 0 |
| `removeNode` | 809 | 809 | 0 | 0 | 0 | 0 |
| `removeByName` | 1014 | 1014 | 0 | 0 | 0 | 0 |
| `get.matrix.format` | 1619 | 1615 | 132 | 4 | 0 | 0 |
| `estimate.nodes.degree` | 1010 | 1010 | 0 | 0 | 0 | 0 |
| `nodes.correlation` | 3039 | 3039 | 480 | 0 | 0 | 0 |
| `nodes.glm` | 829 | 829 | 6 | 0 | 0 | 0 |
| `nodes.fastLm` | 503 | 501 | 3 | 0 | 2 | 0 |
| `SEARCH` | 1218 | 1218 | 382 | 0 | 0 | 0 |
| `write.data.files` | 1221 | 1221 | 210 | 0 | 0 | 0 |
| `.SEARCH.randomise.inner` | 326 | 326 | 92 | 0 | 0 | 0 |
| `SEARCH.randomise` | 492 | 492 | 157 | 0 | 0 | 0 |
| `summarise.iterations` | 1623 | 1623 | 559 | 0 | 0 | 0 |
| **package** | 16526 | 16520 | 2021 | 4 | 2 | 0 |
| `v` | 3 | 3 | 0 | 0 | 0 | 0 |
| `v[names] <- values` | 606 | 605 | 0 | 1 | 0 | 0 |
| `v[names]` | 604 | 604 | 0 | 0 | 0 | 0 |
| `v[v < t]` | 607 | 607 | 0 | 0 | 0 | 0 |
| `v[v > t]` | 607 | 607 | 0 | 0 | 0 | 0 |
| `v[abs(v) > t]` | 607 | 607 | 0 | 0 | 0 | 0 |
| `names(v) <- names` | 608 | 608 | 172 | 0 | 0 | 0 |
| `union(a, b)` | 605 | 605 | 0 | 0 | 0 | 0 |
| `intersect(a, b)` | 605 | 605 | 0 | 0 | 0 | 0 |
| `setdiff(a, b)` | 605 | 605 | 0 | 0 | 0 | 0 |
| `a %in% b` | 605 | 605 | 0 | 0 | 0 | 0 |
| `unique(a)` | 602 | 602 | 0 | 0 | 0 | 0 |
| `c(a, b)` | 605 | 605 | 0 | 0 | 0 | 0 |
| `sort(v, decreasing = TRUE)` | 605 | 605 | 0 | 0 | 0 | 0 |
| `as.numeric(sprintf("%.3f", v))` | 5 | 5 | 0 | 0 | 0 | 0 |
| `a:b` | 7 | 7 | 1 | 0 | 0 | 0 |
| `write.table(x)` | 621 | 621 | 1 | 0 | 0 | 0 |
| `read.table(file)` | 617 | 617 | 5 | 0 | 0 | 0 |
| **rules of R** | 9124 | 9123 | 179 | 1 | 0 | 0 |
