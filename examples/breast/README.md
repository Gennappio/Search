# Breast example (search2014 preset)

`cohorts.yaml` runs the `search2014` preset on Breastdata1 (1231 probesets x 77 samples), with
Breastdata2 loaded as a scoring cohort; the preset has combination and scoring switched off, so
Breastdata2 only goes through the quality checks.

`seeds_breast_cluster.csv` holds a hub probeset and its five strongest neighbours in Breastdata1
(|Spearman| > 0.5), which is the kind of co-expressed seed set the growth loop needs. The five
probesets of the "HIF-1-alpha transcription factor network" of the bundled pathway file are not
co-expressed in these data (|Spearman| <= 0.3 between any two), so with them every seed is pruned
at the first iteration and the run stops with an explicit error; the growth tables under
`results/` show what was pruned.
