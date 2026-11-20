# Reference material (read-only)

Original R sources that Search ports and generalises. They are kept for attribution under the
GPL-2 licence and for consultation; nothing here is executed or imported by Search or its tests.
The only code that runs them is the comparison harness under `validation/r_parity`, which reads
these files as they are (see `docs/porting.md`); a change to a file here changes what the port
is compared with.

- `r_sources/SEARCH/`: the SEARCH R package 0.0.2 (Syed Haider, Francesca M. Buffa, 2014).
- `r_sources/gbm/`: the 2026 R script `run_buffa2010_hpc.R` (two variants) and its SLURM files,
  a single-pass implementation of Buffa et al. 2010 on glioblastoma RNA-seq cohorts.

Known defects of the originals and how Search resolves them are listed in `docs/decisions.md`
and, with the evidence for each, in `docs/porting.md`.
