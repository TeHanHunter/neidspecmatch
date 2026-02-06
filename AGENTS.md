# Repository Guidelines

## Project Structure & Module Organization
Core algorithms sit in `neidspecmatch/` (spectral fitting kernels, priors, configs, and helpers). Scripts such as `run_neidspecmatch.py`, `run_crossval.py`, and `run_neidsm_surfsup.py` demonstrate canonical workflows and make good templates for new automation. Downloaded libraries (CSV metadata, FITS templates, cross-validation caches) belong inside `library/` and are located through `neidspecmatch.config`. `tests/` stores reference plots and NumPy vectors for regression checks, while notebooks under `tutorial/` provide curated walkthroughs.

## Build, Test, and Development Commands
Editable install for day-to-day work:
```bash
pip install -e .
```
Publishable artifacts:
```bash
python -m build
```
Smoke test on an example FITS (library must be populated):
```bash
python run_neidspecmatch.py
```
Recompute cross-validation statistics for all orders:
```bash
python run_crossval.py
```

## Coding Style & Naming Conventions
Default to PEP 8 with 4-space indentation. Use snake_case for functions and variables, ALL_CAPS for module constants, and CamelCase only for classes (see `target.py`). Group imports stdlib/third-party/local, and prefer keyword arguments when calling helpers with many parameters (e.g., `run_specmatch_for_orders`). You may run `ruff format` or `black` locally, but avoid repo-wide style-only commits. Store shared settings in `config.py` and expose runtime knobs through explicit parameters rather than hidden globals.

## Testing Guidelines
We do not have CI yet, so add deterministic pytest modules in `tests/` whenever you touch model logic. Use the `test_<feature>.py` pattern and light fixtures (CSV snippets, trimmed FITS, PNG baselines) saved beside the tests. Record the Teff/[Fe/H]/log g scatter you observe and stash plots in `tests/img`. Always run the relevant workflow script—`python run_crossval.py` for library or statistics edits, `python run_neidspecmatch.py` for fitting changes—and include the headline metrics in your PR.

## Commit & Pull Request Guidelines
Recent commits mix concise imperatives (`add figure`) with one-line summaries (`Refactor NEID spectrum processing…`). Keep using that style: subject under 72 characters, present-tense imperative, with optional wrapped detail lines. Pull requests must explain motivation, enumerate affected orders or modes, link related Zenodo/Drive files, and state which commands from this guide were run. Attach or link updated diagnostic plots (`tests/img`) whenever a visualization or uncertainty estimate shifts.

## Security & Configuration Tips
Never commit Zenodo payloads or proprietary spectra—document the download links and depend on the paths in `config.py`. If a workflow needs credentials or remote paths, read them from environment variables with safe fallbacks. Validate user-supplied FITS paths before iterating to avoid directory traversal, and redact observer metadata from logs or plots before sharing externally.***
