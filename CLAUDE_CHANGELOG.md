# Claude Changelog

Record of all changes made by Claude Code in this repository.

---

## 2026-04-23

- **Added GitHub Actions workflow** (`.github/workflows/test.yml`) — runs the test suite on every push and pull request to `devel` and `main`. Uses the official `bioconductor/bioconductor_docker:RELEASE_3_22` container. Caches R packages (keyed to `DESCRIPTION`) and ExperimentHub data (static key since datasets EH7858–EH7863 are versioned) to minimise run time.

- **Added validation scripts** (`validation/`) — three scripts for comparing Bioconductor release vs dev outputs to verify correctness of the CV_Tissue merge rewrite: `01_run_bioc.R`, `02_run_dev.R`, `03_compare.R`. Tests single-species (must match exactly), two-species, and three-species cases. Includes manual spot-check instructions.
- **Updated `.Rbuildignore`** — added `validation/`, `CLAUDE.md`, and `CLAUDE_CHANGELOG.md` so they are excluded from the package tarball.

## 2026-04-22

- **Created `CLAUDE.md`** — initial guidance file documenting project overview, common commands, architecture, key files, and conventions.
- **Removed deprecation markers** — CoSIA was listed on a Bioconductor deprecation candidate list but was never actually removed; it remains active in Bioconductor 3.22 at v1.10.1. Removed the `.Deprecated()` call from `R/zzz.R` and `PackageStatus: Deprecated` from `DESCRIPTION`.
- **Created unit test suite** — added `tests/testthat.R` and five test files covering the full public API:
  - `tests/testthat/helper-CoSIA.R` — shared fixtures and `skip_if_no_experimenthub()` helper
  - `tests/testthat/test-CoSIAn-constructor.R` — 11 tests for constructor validation (offline)
  - `tests/testthat/test-viewCoSIAn.R` — 4 tests for the `viewCoSIAn` accessor (offline)
  - `tests/testthat/test-getConversion.R` — 5 tests for gene ID/ortholog conversion; 4 run offline via annotationDBI + HomoloGene, 1 skips without network (NCBIOrtho)
  - `tests/testthat/test-getGEx.R` — 5 tests for ExperimentHub expression data retrieval (skip if offline)
  - `tests/testthat/test-getGExMetrics.R` — 9 tests covering all 6 metric types (skip if offline)
- **Updated `DESCRIPTION`** — added `Config/testthat/edition: 3` to declare testthat edition 3.
- **Updated `CLAUDE.md`** — added Claude Changelog section requiring this file to be updated with every Claude session; corrected stale note about `zzz.R`.
- **Fixed `viewCoSIAn` bug** (`R/CoSIA-utility-method.R`) — replaced `ifelse()` with `if/else`. `ifelse()` is vectorized and strips S3 attributes, causing it to return a bare scalar instead of the data frame slot value.
- **Fixed `CV_Tissue` multi-species bug** (`R/CoSIA-getGExMetrics-method.R`) — rewrote the merge loop that accumulated the growing result back into itself on each iteration, creating a Cartesian product with cascading `.x`/`.y` duplicate columns and a broken `Anatomical_entity_name` column (also contained a typo `"Anatomical_enti_name.y"`). Replaced with a per-species independent left join against the original `id_dataframe`, combining results by joining on all ID columns + `Anatomical_entity_name`.
