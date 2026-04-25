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

# PESTO 0.1.0.9000 (development version)

## New vignettes

* `pestpp-comparison-and-simulation` — three-part vignette covering
  (i) low-dimensional well-posed comparison of PESTO native IES versus
  the upstream `pestpp-ies` binary on a shared analytical problem,
  (ii) a 100-parameter Tikhonov-regularised inverse problem with SVD
  truncation sensitivity and rank-dependent rSVD/LAPACK benchmarking,
  and (iii) a 50-replicate Monte-Carlo simulation study that exercises
  every exported function in `NAMESPACE`.

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
