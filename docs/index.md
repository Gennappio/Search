# Search

Seed-based gene co-expression clustering with named presets: one engine, the `buffa2010` and
`search2014` presets, bootstrap and null models, cohort combination, signature scores.

- [Method](methods.md): the formulas of every stage.
- [Presets](presets.md): what each preset reproduces.
- [Decision log](decisions.md): deviations from the R originals and the open questions.
- [Command line](cli.md): commands, the cohort manifest and the run directory.
- [Validation](validation.md): the tests of every layer and the results on real data.

## Quick start

```bash
uv sync --all-extras
uv run search run --manifest examples/gbm/cohorts.yaml --preset buffa2010 --out runs --n-jobs 8
uv run search plot runs/latest
```

`examples/gbm/cohorts.yaml` expects the GBM data folder next to it; edit the paths or copy the
file next to the data. `examples/breast/cohorts.yaml` runs the `search2014` preset on the breast
cohorts bundled as test fixtures. The documentation site is not published (GitHub Pages is not
available for private repositories on the free plan): read it with `uv run mkdocs serve`. The optional pysigscore scorers need the package installed in the
environment by hand (it is not on PyPI): `uv pip install --python .venv/bin/python <path-or-git-url>`.
