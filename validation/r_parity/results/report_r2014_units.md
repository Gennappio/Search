Unit cases of `search.r2014` against R version 4.4.3 (2025-02-28) on aarch64-apple-darwin20, 2026-11-20.

| function or rule of R | cases | equal | of which R stops | refused by the port, as stated | differ | skipped |
|---|---|---|---|---|---|---|
| `getO` | 73 | 73 | 0 | 0 | 0 | 0 |
| `getOmax` | 90 | 90 | 0 | 0 | 0 | 0 |
| `removeNode` | 49 | 49 | 0 | 0 | 0 | 0 |
| `removeByName` | 64 | 64 | 0 | 0 | 0 | 0 |
| `get.matrix.format` | 99 | 95 | 11 | 4 | 0 | 0 |
| `estimate.nodes.degree` | 60 | 60 | 0 | 0 | 0 | 0 |
| `nodes.correlation` | 189 | 189 | 26 | 0 | 0 | 0 |
| `nodes.glm` | 70 | 70 | 5 | 0 | 0 | 0 |
| `nodes.fastLm` | 48 | 48 | 4 | 0 | 0 | 0 |
| `SEARCH` | 78 | 78 | 29 | 0 | 0 | 0 |
| `write.data.files` | 81 | 81 | 18 | 0 | 0 | 0 |
| `.SEARCH.randomise.inner` | 22 | 22 | 7 | 0 | 0 | 0 |
| `SEARCH.randomise` | 36 | 36 | 13 | 0 | 0 | 0 |
| `summarise.iterations` | 103 | 103 | 40 | 0 | 0 | 0 |
| **package** | 1062 | 1058 | 153 | 4 | 0 | 0 |
| `v` | 3 | 3 | 0 | 0 | 0 | 0 |
| `v[names] <- values` | 36 | 35 | 0 | 1 | 0 | 0 |
| `v[names]` | 34 | 34 | 0 | 0 | 0 | 0 |
| `v[v < t]` | 37 | 37 | 0 | 0 | 0 | 0 |
| `v[v > t]` | 37 | 37 | 0 | 0 | 0 | 0 |
| `v[abs(v) > t]` | 37 | 37 | 0 | 0 | 0 | 0 |
| `names(v) <- names` | 38 | 38 | 16 | 0 | 0 | 0 |
| `union(a, b)` | 35 | 35 | 0 | 0 | 0 | 0 |
| `intersect(a, b)` | 35 | 35 | 0 | 0 | 0 | 0 |
| `setdiff(a, b)` | 35 | 35 | 0 | 0 | 0 | 0 |
| `a %in% b` | 35 | 35 | 0 | 0 | 0 | 0 |
| `unique(a)` | 32 | 32 | 0 | 0 | 0 | 0 |
| `c(a, b)` | 35 | 35 | 0 | 0 | 0 | 0 |
| `sort(v, decreasing = TRUE)` | 35 | 35 | 0 | 0 | 0 | 0 |
| `as.numeric(sprintf("%.3f", v))` | 5 | 5 | 0 | 0 | 0 | 0 |
| `a:b` | 7 | 7 | 1 | 0 | 0 | 0 |
| `write.table(x)` | 51 | 51 | 1 | 0 | 0 | 0 |
| `read.table(file)` | 47 | 47 | 5 | 0 | 0 | 0 |
| **rules of R** | 574 | 573 | 23 | 1 | 0 | 0 |
