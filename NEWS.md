# PESTO 0.4.1

## AAGI guidance uplift

A code-aesthetics and review-readability patch on top of `0.4.0`. No
runtime behaviour changes; no exported-API changes; no shipped-data
changes. The aim is to lift every source surface to the bar set out in
`r_style.md` ahead of any AAGI-AUS push.

### Validation delegation (r_style.md invariant 9)

* `R/internal_validation.R` introduces the shared primitive validators
  (`.assert_positive_scalar()`, `.assert_nonneg_scalar()`,
  `.assert_character_scalar()`, `.assert_logical_scalar()`,
  `.assert_path_exists()`, `.assert_matrix()`,
  `.assert_numeric_vector()`, `.assert_function()`,
  `.assert_data_frame()`, `.assert_choice()`, `.assert_same_ncol()`,
  `.assert_same_nrow()`, `.assert_required_cols()`). All
  `@noRd @keywords internal`; every helper signals failure via
  `stop(call. = FALSE, ...)` with a backticked argument name.
* Public functions in `apsim_callback.R`, `pesto_reference_ies.R`,
  `pesto_run.R` (`pesto_ies_callback`), `pst_io.R`, `scenario.R`,
  `surrogate.R`, `manifest.R`, `ensemble_io.R`, `plot.R`, and
  `check_surrogate_regime.R` now open with `.check_*` / `.assert_*`
  calls instead of inline `if (!is.x) stop(...)` walls.
* Error messages backtick the offending argument name throughout
  (Sparks convention).

### Section banners (r_style.md invariant 4)

* Long function bodies (`pesto_ies_callback`, `pesto_ies`,
  `pesto_glm`, `pesto_sweep`, `pesto_sensitivity`, `read_pst`,
  `write_pst`, `apsim_callback`, `pesto_reference_ies`,
  `.find_pestpp_exe`) now carry Sparks-style dash-banner section
  comments that paragraph the work (validate inputs / resolve paths /
  iterate / parse outputs / assemble result).

### Import concentration (r_style.md invariant 10)

* `@importFrom ggplot2` annotations consolidated into
  `R/pesto-package.R`; the per-function `@importFrom` annotation on
  `plot_phi()` has been removed in favour of inline `ggplot2::`
  qualification at the call sites.

### Vignette prose

* Prose semicolons in `vignettes/apsim-callback.Rmd` and
  `vignettes/ensemble-manifest.Rmd` converted to `. Capital` joins per
  `manuscript_style.md` invariant 5.

### README

* Added `Dependencies`, `Contributing`, and `Acknowledgements`
  sections per the AAGI repository-guidelines README contract.
  Citation block bumped to `R package version 0.4.1`.

# PESTO 0.4.0

## AAGI recipes uplift and canon channel migration

This release contains no R, C++, or shipped-data changes. It is a
governance, metadata, and project-hygiene release that lands the AAGI
canon recipes on the `max578/PESTO` channel.

### Canon channel migration

* Primary canon channel migrated from `AAGI-AUS/PESTO` to
  `max578/PESTO`. `DESCRIPTION`, `CITATION.cff`, `codemeta.json`,
  `inst/CITATION`, `_pkgdown.yml`, `README.md`, `CONTRIBUTING.md`,
  `API_STABILITY.md`, and the pkgdown GitHub Actions workflow header
  now point to `https://github.com/max578/PESTO` and
  `https://max578.github.io/PESTO`. The `aagi` git remote is retained
  as a frozen read-only mirror; no push to `AAGI-AUS` without explicit
  per-instance maintainer approval.
* Package-root `CLAUDE.md` declares `aagi_aus: out-of-scope` so the
  AAGI-AUS canon signal-detection deactivates for this package. The
  file is excluded from R-package builds via `.Rbuildignore`.
* `man/PESTO-package.Rd` regenerated to inherit the new URLs from
  `DESCRIPTION` via `devtools::document()`.

### Metadata version sync

* `CITATION.cff` (`version: 0.1.0`), `codemeta.json` (`version: "0.1.0"`),
  `inst/CITATION` (`R package version 0.3.3`), and the README citation
  block were not in lock-step with `DESCRIPTION`. All four are now on
  `0.4.0` with `date-released: "2026-05-28"` and
  `dateModified: "2026-05-28"`.

### Sole copyright holder

* `codemeta.json` `copyrightHolder` corrected from
  `Organization "Supremum Consulting Ltd"` to
  `Person "Max Moldovan"` with ORCID `0000-0001-9680-8474` and
  University of Adelaide affiliation, matching `Authors@R` and
  `LICENSE.md`.

### AAGI canon recipes added

* `CODE_OF_CONDUCT.md` (Contributor Covenant v2.1, pointer form).
* `SECURITY.md` (vulnerability reporting policy; maintainer email,
  five-working-day acknowledgement, scope statement).
* `air.toml` (Air formatter configuration: 80-char line width,
  two-space indent, auto line endings).
* `.lintr` (lintr defaults aligned with `r_style.md` direction:
  80-char line, `snake_case` / `dotted.case` / symbols object names,
  two-space indent; `src`, `tools`, `inst/extdata`, `vignettes`
  excluded).
* `.Rbuildignore` already excluded all four paths; no tarball impact.

# PESTO 0.3.3

## FLIBS Makevars portability fix (closes critical-review P2 #7)

* `src/Makevars` `PKG_LIBS` now follows `$(BLAS_LIBS)` with
  `$(FLIBS)` per Writing R Extensions §1.2.1.5. Resolves the
  pre-existing structural WARNING
  ("apparently using $(BLAS_LIBS) without following $(FLIBS) in
  'src/Makevars'") that had survived several check passes as a
  "documented local-env artifact". `Makevars.win` was already
  correct; no change needed there.
* The remaining baseline structural WARNING (Apple-clang 21.0.0 vs.
  `R_ext/Boolean.h:62` `-Wfixed-enum-extension` pragma) is in R's
  own header on this toolchain version and is not present on CRAN's
  build farm; it persists harmlessly.

# PESTO 0.3.2

## CSV-only manifest mode renamed (post critical-review, 2026-05-16)

* `write_manifest(format = "csv")` is renamed to
  `write_manifest(format = "csv_unverified")` to flag the weaker
  integrity contract at every call-site. The mode itself is unchanged
  (CSV-only sidecars, hash recorded but not disk-verifiable).
* New top-level YAML field `integrity: verifiable | not_verifiable`
  derived from `format`. Verifiable: `rds`, `both`. Not verifiable:
  `csv_unverified`. Lets non-R downstream tools (paper-skill graders,
  Python pipelines) branch on the integrity contract without parsing
  the PESTO-specific format vocabulary.
* The legacy `format = "csv"` spelling is still accepted at the API
  boundary with a deprecation warning; the persisted form always uses
  `csv_unverified`. `read_manifest()` normalises old YAMLs on the
  read side, so 0.3.1 manifests round-trip cleanly under 0.3.2.
* Validator vocabulary updated:
  `{rds, both, csv_unverified}`. The old `"csv"` token is rejected at
  slot-set time (only the renamed argument accepts it, with a warning).
* New tests cover the renamed value, the deprecation-warning path,
  and the `integrity:` YAML field for verifiable modes.
* Vignette `ensemble-manifest.Rmd` reframed: explicit "Inspection
  CSVs (verifiable, via format = 'both')" vs. "Unverified CSV export
  (via format = 'csv_unverified')" sections; the latter is presented
  as "for export, not for storage you intend to re-load and trust".

# PESTO 0.3.1

## Manifest sidecar `format=` option (roadmap §A5 polish)

* `write_manifest()` gains a `format = c("rds", "both", "csv")`
  argument. `"rds"` (default) preserves the current bit-exact binary
  behaviour. `"both"` writes RDS sidecars plus parallel CSV inspection
  files (`*_inspection.csv`); the SHA-256 hash stays bound to the RDS.
  `"csv"` writes CSV-only sidecars for inspection / interchange
  workflows where bit-exact integrity is not required.
* New S7 slot `format` on `pesto_ensemble_manifest` records the
  on-disk serialisation mode (default `"rds"`; preserved through
  read/write round-trips). Validator enforces the three-value vocabulary.
* `read_manifest()` dispatches on file extension in the YAML's
  `artefacts:` block — reads RDS via `readRDS()`, CSV via
  `utils::read.csv()`.
* `verify_manifest()` gains a `message` field on its return list and
  returns `ok = NA` (with explanation) for `format = "csv"` manifests
  whose IEEE 754 doubles have round-tripped through a write formatter.
  Existing `format = "rds"` callers see no behaviour change; the new
  field is `NULL` in that case.
* YAML schema picks up a top-level `format:` key and an optional
  `inspection_csv:` block when `format = "both"`. Backwards-compatible:
  YAMLs written by PESTO 0.3.0 (no `format:` key) read back with
  `format = "rds"` per the default. No schema-version bump required.
* New tests in `test-manifest.R` cover all three formats plus the
  unknown-format rejection path.

# PESTO 0.3.0

## Ensemble Manifest as S7 Cross-Package Contract (roadmap §A5)

* New S7 class `pesto_ensemble_manifest` — versioned, hashed,
  provenance-tracked container for ensemble-run output. Slots cover
  `params`, `outputs`, `weights`, `obs_target`, `seed`, `data_hash`
  (SHA-256), `fidelity`, `apsim_version`, `pesto_version`, `timestamp`,
  plus method context (`method`, `noptmax`, `lambda_schedule`,
  `failure_rate`). This is the contract object that downstream
  consumers (`kernR`, `proxymix`, paper-skill) will read.
* New `as_manifest()` — S7 generic with a method for
  `pesto_ies_callback_result`. Non-destructive: wraps without mutating
  the source result.
* New `write_manifest()` / `read_manifest()` — YAML+RDS serialisation.
  The YAML carries metadata + relative paths to three sidecar RDS files
  (`*_params.rds`, `*_outputs.rds`, `*_assim.rds`); RDS is used in
  preference to CSV so IEEE 754 doubles round-trip bit-exactly (the
  SHA-256 integrity check would otherwise trip on CSV formatter
  precision loss).
* New `verify_manifest()` — recomputes the SHA-256 over
  `(params, outputs, weights, obs_target, seed)` and compares to the
  stored value, returning a diagnostic list. Detects post-write
  tampering with the sidecar CSVs.
* New vignette: `ensemble-manifest` — end-to-end demo of construct →
  write → read → verify, plus tamper-detection.

## In-Process IES via R Callback (UQ ag-stack roadmap §D4)

* New `pesto_ies_callback()` — drives an Iterative Ensemble Smoother
  entirely in R using a user-supplied forward-model callable, bypassing
  the `.pst`-file write/read cycle of `pesto_ies()`. Each iteration calls
  the existing C++ kernel `ensemble_solution()` (Chen & Oliver, 2013).
  Tolerates per-realisation failures via `on_failure = c("na", "stop")`
  and reports a `failure_rate` in the result object. Phase-1 behaviour
  uses a single lambda per iteration (or user-supplied schedule); a
  full `pestpp-ies`-style lambda line-search is a planned Phase-2
  enhancement.
* New `apsim_callback()` — adapter that wraps the `apsimx` package (now
  in `Suggests`) into a forward-model closure suitable for
  `pesto_ies_callback()`. Per-realisation template copy, parameter edit
  via `apsimx::edit_apsimx()` / `edit_apsim()`, run, and extraction.
  Failures (edit / run / extractor) surface as `NA` rows for the IES
  driver to handle.
* New vignette: `apsim-callback` — synthetic linear-Gaussian recovery
  demo plus disabled apsimx example.
* New benchmark script: `inst/benchmarks/d4_callback_vs_pst.R`.
* `pesto_ies_callback()` records `obs_target`, `obs_sd`, and `weights`
  on its result list so the manifest emitter has full IES context to
  capture.

## Notes

* Real-APSIM benchmarking (target ≥10× speed-up vs `.pst` path) is
  deferred to the §D1 scenario library landing. The current benchmark
  script measures the callback path on a synthetic surrogate forward
  model only.
* `apsimx` is in `Suggests` not `Imports`; `apsim_callback()` checks
  for it at call time with `requireNamespace()`.
* New hard imports: `S7 (>= 0.2.0)`, `yaml (>= 2.3.0)`,
  `digest (>= 0.6.0)`.

# PESTO 0.2.0 (2026-04-25)

## Bug fixes

* **I1 — `ensemble_solution()` sign-convention bug.** The C++ kernel
  requires `obs_resid = sim - obs`; the docstring previously stated
  the inverse. Two genuine in-package call sites were silently
  inverting upgrades: `src/surrogate_ies.cpp:347` and
  `vignettes/surrogate-ies.Rmd:148, 181`. Both fixed; surrogate-IES
  now applies upgrades in the correct direction. Regression test
  `tests/testthat/test-ensemble-solution-sign.R` asserts strict
  monotone phi descent under the correct convention AND geometric
  divergence under the inverted one.

## New exported functions

* `pesto_reference_ies()` — pure-R, textbook implementation of the
  Chen & Oliver (2013) eq. 12 IES update. Independent of the C++
  kernel; used as the canonical comparison target by the comparison
  vignette so it ships and runs without the upstream `pestpp-ies`
  binary. Cross-validated against the C++ kernel at machine precision
  (max element-wise delta = 5.8e-15).
* `check_surrogate_regime()` — soft guardrail that warns when the
  surrogate-IES regime is unfavourable (`n_train < threshold * n_params`).
  Stand-alone helper, not auto-invoked by `pesto_surrogate_ies()`
  (v0.3 wiring candidate).

## API enhancements

* `plot_identifiability()` gains a `jacobian = NULL` matrix-input
  path. Backward-compatible: `jco_file = NULL` retained; the two are
  mutually exclusive.

## Self-contained `pestpp-ies` comparison

* `vignettes/pestpp-comparison-and-simulation.Rmd` now compares PESTO
  native IES against `pesto_reference_ies()` by default — no upstream
  binary required. The pure-R reference cache ships at
  `inst/extdata/pestpp_cache/scenario_a_reference.rds` (SHA-256-pinned
  to the prior ensemble). When the developer-side cache
  `tools/pestpp_benchmark/scenario_a_pestpp_ies.rds` is present and
  `PESTO_PESTPP_BIN` resolves, the vignette extends the agreement
  plot with the live binary's posterior.
* Hardcoded `/Users/a1222812/...` path replaced with
  `Sys.getenv("PESTO_PESTPP_BIN")` + `Sys.which("pestpp-ies")`
  fallback.
* New `tools/pestpp_benchmark/run_benchmark.R` regenerates both
  caches deterministically. Documented in `CONTRIBUTING.md`.

## Documentation and CRAN-readiness

* Every export now has a runnable `@examples` block (30 of 30
  documented exports/methods). The four external-binary runners
  (`pesto_ies`, `pesto_glm`, `pesto_sweep`, `pesto_sensitivity`) and
  `pesto_surrogate_ies` use guarded `\donttest{}` (no `\dontrun{}`).
* Vignettes acquire a "Regime of applicability" subsection
  (`surrogate-ies.Rmd`) and an "Honest reading — surrogate savings in
  this regime" defence paragraph (`pestpp-comparison-and-simulation.Rmd`
  Section 3) covering the curse-of-dimensionality finding from
  investigation I3.
* Kernel docstring `?ensemble_solution` now states the `sim - obs`
  convention with a full GLM-derivation rationale.
* New `cran-comments.md` with per-NOTE justification.
* New `CITATION.cff` (CFF 1.2.0 + ORCID + preferred citation).
* New `codemeta.json` (CodeMeta 2.0).
* New `CONTRIBUTING.md` documenting the developer benchmark workflow.
* New `inst/WORDLIST` with ~120 domain terms; `Language: en-AU`
  added to DESCRIPTION.

## Build and CRAN-portability

* `src/Makevars`: `PKG_LIBS` gains `$(FLIBS)` (CRAN portability
  requirement for `$(BLAS_LIBS)`).
* `LICENSE` renamed to `LICENSE.md` and `.Rbuildignore`-d (CRAN
  convention for GPL-licensed packages).

## Open investigations (status)

* I1 — **CLOSED** (this release).
* I2 — partially closed: PESTO C++ kernel verified equivalent to
  Chen & Oliver (2013) eq. 12 at machine precision. Residual
  PESTO-vs-pestpp-ies posterior gap is structural (Marquardt
  sub-cycling + per-realisation noise perturbation in upstream
  pestpp-ies) and is a v0.3 implementation milestone.
* I3 — **DOCUMENTED** (this release): regime-of-applicability text
  added to both vignettes; `check_surrogate_regime()` helper exported.

# PESTO 0.1.1

## Sole-Authorship Consolidation

* Authorship consolidated to Max Moldovan as sole `aut`, `cre`, and `cph`.
  Supremum Consulting Ltd. removed from `Authors@R` (administrative
  consolidation by sole director; no licence change).
* Licence unchanged: GPL-3 or any later version.
* `LICENSE` file rewritten with corrected canonical wording and copyright
  attribution.
* Source-file headers in `src/` updated to reflect sole authorship.

# PESTO 0.1.0

## Initial Release

### Core Features
* `ensemble_solution()` — High-performance C++ implementation of the IES
  ensemble update equation (Chen & Oliver, 2013) via RcppEigen.
* `ensemble_solution_mda()` — Multiple Data Assimilation (Evensen, 2018)
  update kernel.
* `compute_phi()` — Fast weighted sum-of-squares objective function.

### Hardware Acceleration
* `adaptive_svd()` — Automatic SVD backend selection (LAPACK, Eigen BDCSVD,
  or randomised SVD) based on matrix size and target rank.
* `rsvd()` — Randomised SVD (Halko-Martinsson-Tropp, 2011) for asymptotically
  faster rank-k approximations.
* `accelerate_svd()` — Direct LAPACK SVD leveraging platform-optimised BLAS
  (Apple Accelerate/AMX on macOS, MKL or OpenBLAS on Linux).
* `ensemble_solution_gpu()` — GPU-ready ensemble solution with adaptive SVD
  backend and performance diagnostics.

### Surrogate-Accelerated IES (Novel)
* `train_gp_surrogate()` — Gaussian Process surrogate model training with
  automatic hyperparameter selection (median heuristic).
* `predict_gp_surrogate()` — GP prediction with uncertainty quantification.
* `surrogate_ensemble_update()` — Surrogate-accelerated IES update with
  adaptive model/surrogate switching and control-variate bias correction.
* `adaptive_ensemble_size()` — Convergence-aware dynamic ensemble sizing
  based on ESS and coefficient of variation diagnostics.

### PEST++ Integration
* `read_pst()` / `write_pst()` — PEST control file I/O.
* `read_ensemble()` / `write_ensemble()` — Ensemble file I/O (CSV + binary).
* `pesto_ies()`, `pesto_glm()`, `pesto_sweep()`, `pesto_sensitivity()` —
  High-level wrappers for PEST++ executables.
* `create_pest_scenario()` — Programmatic scenario builder.

### Visualisation
* `plot_phi()` — Objective function convergence plotting.
* `plot_ensemble()` — Prior/posterior parameter distribution comparison.
* `plot_identifiability()` — SVD-based parameter identifiability analysis.
* `plot_surrogate_diagnostics()` — Surrogate IES performance visualisation.
