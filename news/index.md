# Changelog

## PESTO 0.10.1

Follows 0.10.0’s PEST++ rebuild by applying the same question to the
other external authority PESTO drives – APSIM – and finding the same
shape of defect: a working engine PESTO could not see, and an adapter
whose tests only ever asked PESTO about itself.

### Bug fixes

- PESTO now finds an APSIM Next Gen engine via `APSIM_EXE_PATH`. It
  consulted only `apsimx`’s `exe.path` option, so a machine with a
  working APSIM and the environment variable exported still reported no
  engine, and
  [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  stamped the closure’s `apsim_version` as `NA` – leaving every run it
  fed unattributable to the simulator that produced it. APSIM has no
  discoverable install location on macOS:
  [`apsimx::apsim_version()`](https://rdrr.io/pkg/apsimx/man/apsim_version.html)
  derives the version from `/Applications` folder names, so it sees
  neither a `~/Applications` install nor a source build, and never reads
  `exe.path` at all. An explicit `apsimx_options(exe.path=)` still takes
  precedence. The adapter now has tests that run against a real engine
  rather than stubs.

## PESTO 0.10.0

This release rebuilds the PEST++ invocation layer. PESTO’s shell-out to
PEST++ had never run: control variables were passed as `/h :name=value`
command-line switches, which PEST++ has never accepted, so every call
exited before the binary started. The defects below have been present
since the initial release (9 April 2026) and shipped in every release
since; they surfaced now because this is the first time the layer was
checked against the PEST++ sources and then against the real binary,
rather than against PESTO’s own tests.

Users of
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
– PESTO’s native R ensemble smoother, and the path the benchmark suite
exercises – are unaffected: it does not shell out.

### Breaking changes

- [`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md)
  and
  [`pesto_glm()`](https://max578.github.io/PESTO/reference/pesto_glm.md)
  now default `noptmax` to `NULL`, meaning “leave the control file’s own
  iteration cap alone”. The previous defaults (`4` and `20`) were
  documented as overrides but, because the whole command-line path was
  broken, had never once been applied to a real run. Honouring them now
  would have silently retuned every caller’s control file, so the safer
  reading of a value the user never chose is to respect the file. Pass
  `noptmax` explicitly to override. A call that injects nothing writes
  no file at all and runs the caller’s own control file, rather than a
  `<base>_pesto.pst` copy of it.

### Bug fixes

- PESTO now finds a PEST++ install that is not on the `PATH`. It
  consulted only the `PATH` and a copy bundled at `inst/bin` – which
  PESTO has never shipped, so that branch could not fire and the
  documented “uses the bundled binary” fallback did not exist. PEST++
  has no installer and is not on the `PATH` by default, so a machine
  could have it installed and configured and still be told it was
  missing. Resolution is now `exe`, then the per-tool environment
  variable (e.g. `PESTPP_IES_EXE_PATH`), then `PESTPP_BIN_DIR`, then the
  `PATH`.

- [`pesto_version()`](https://max578.github.io/PESTO/reference/pesto_version.md)
  reports a version rather than a failed run’s log, and no longer writes
  to the caller’s working directory. PEST++ has no `--version` flag: it
  read `--version` as a control-file name, printed its banner, failed on
  the missing `--version.pst`, and left `--version.log`, `--version.rec`
  and `--version.rst` behind. `$pestpp_version` held that whole
  transcript, `std::exception` and all, and looked plausible only
  because the banner is printed before the error. The banner is now
  parsed for its `version:` line, from a temporary directory.

- The PEST++ invocation layer has been rebuilt. PESTO passed control
  variables as `/h :name=value` command-line switches, which PEST++ has
  never accepted: its parser takes a control file, an optional `/r` or
  `/j`, and an optional run-manager switch, where `/h` selects the
  PANTHER run manager and expects a `host:port`. Every
  [`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md)
  call built at least one such switch, so the binary exited with a
  command-line error before starting. Control variables now go where
  PEST++ reads them: `noptmax` into `* control data`, and every other
  option as a `++key(value)` line. No test caught this because the
  binary was never a test dependency; the suite checked PESTO against
  PESTO. Verified against the PEST++ sources rather than against PESTO’s
  own belief about them.

- [`create_pest_scenario()`](https://max578.github.io/PESTO/reference/create_pest_scenario.md)
  now writes instruction files into the control file.
  `instruction_files` was read only to count them, so every control file
  PESTO produced declared `NINSFLE` and then supplied none, and PEST++
  refused all of them:
  `model input/output error: number of instruction files = 0`. PESTO had
  never written a runnable control file; nothing local could see it,
  because the control file is only ever read by the binary.

- `$exit_code` is an exit code again. With `verbose = FALSE` – which
  every example uses –
  [`system2()`](https://rdrr.io/r/base/system2.html) returns the
  captured OUTPUT rather than the status, so `$exit_code` held PEST++’s
  entire log.

- [`pesto_glm()`](https://max578.github.io/PESTO/reference/pesto_glm.md)
  now honours `noptmax` and `extra_args`, and
  [`pesto_sensitivity()`](https://max578.github.io/PESTO/reference/pesto_sensitivity.md)
  now honours `extra_args`. All three were documented, exported, and
  read by nothing.

- `pesto_sensitivity(method = "sobol")` now runs Sobol. `method`
  selected the label on the returned object but never reached the
  binary, so pestpp-sen ran its Morris default and the result was
  reported as Sobol regardless.

- PEST++ now runs in the directory holding the control file. The file is
  passed by basename, so it – and every relative template, instruction,
  and model-command path inside it – was previously resolved against R’s
  working directory, which is only correct by coincidence.

- [`read_pst()`](https://max578.github.io/PESTO/reference/read_pst.md)
  now reads `NOPTMAX`, and
  [`write_pst()`](https://max578.github.io/PESTO/reference/write_pst.md)
  emits it instead of a hard-coded `30`. A read/write round-trip
  previously reset the caller’s iteration cap without saying so.

- [`write_ensemble()`](https://max578.github.io/PESTO/reference/write_ensemble.md)
  rejects a `format` it cannot write. `format = "binary"` wrote a CSV
  under the requested name.

- [`plot_identifiability()`](https://max578.github.io/PESTO/reference/plot_identifiability.md)
  honours `pst`, documented as the source of parameter names for a
  `.jco` that carries none but never consulted. Reading such a file also
  dropped columns: blank labels all collided on one name, so a
  three-column Jacobian came back with one.

- Runs now write `<name>_pesto.pst` beside the control file and name
  their outputs after it. This is the exact input PEST++ was given,
  recorded for reproducibility.

### Minor improvements and fixes

- [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  now reads the engine version for the `"apsim_version"` provenance
  attribute from `Models --version` on the configured binary, rather
  than from
  [`apsimx::apsim_version()`](https://rdrr.io/pkg/apsimx/man/apsim_version.html).
  On macOS the latter infers the version from `/Applications` folder
  names and cannot see a source build or an install under
  `~/Applications`, so the attribute previously fell back to `NA`; it
  now reflects the simulator actually in use.

## PESTO 0.9.0

- **Vignettes brought to publication quality.** *Benchmarking PESTO
  against PEST and PEST++* (formerly the “comparison and simulation
  study”) is reduced from roughly 1300 lines to a focused benchmark –
  comparative estimates, the notable deviations explained, and compute
  times – plus one worked example; its coverage-driving simulation study
  moved into `tests/` (`test-export-surface.R`). The pkgdown site now
  renders maths with KaTeX (resolving a `\boldsymbol` extension-load
  failure), all vignettes use sentence-case headings and verified
  references (the PEST citation is corrected to the 2015 *Calibration
  and Uncertainty Analysis* book), and
  [`plot_identifiability()`](https://max578.github.io/PESTO/reference/plot_identifiability.md)
  gains a `top_n` cap with a ranked-lollipop layout so high-dimensional
  problems stay legible. A `tools/check_publication_quality.R` guard
  blocks internal-flag words from shipped prose.
- **Name and title grounded in the lineage.** PESTO expands as
  *Parameter ESTimation Optimised* – the **PEST** approach (PEST =
  *Parameter ESTimation*, Doherty 2015) brought to R and optimised –
  rather than the earlier letter-by-letter backronym. The package title
  is now *Parameter Estimation Optimised, with APSIM Coupling*, and the
  website tagline, citation, and GitHub description match.
- **A “Performance characteristics” section** in *Benchmarking PESTO
  against PEST and PEST++* documents the four levers behind the speed
  and evaluation economy: the in-process cost model (wall-clock speed-up
  at matched accuracy), surrogate evaluation-economy, the multi-fidelity
  control-variate variance reduction (Kennedy & O’Hagan 2000; Glasserman
  2003), and convergence-based early stopping – each with its maths and
  a measured outcome, and with the limits stated plainly (not fewer
  solves than classic PEST; raw intervals under-cover until inflation).
  The README’s surrogate feature is corrected from a flat “50-90%” to
  the measured, regime-dependent saving.
- **[`ensemble_solution_gpu()`](https://max578.github.io/PESTO/reference/ensemble_solution_gpu.md)
  renamed to
  [`ensemble_solution_adaptive()`](https://max578.github.io/PESTO/reference/ensemble_solution_adaptive.md).**
  The former name implied GPU computation; the function performs CPU
  adaptive-SVD backend selection (randomised SVD versus a dense LAPACK /
  Accelerate decomposition), and its documentation no longer claims CUDA
  / cuSOLVER support.
  [`ensemble_solution_gpu()`](https://max578.github.io/PESTO/reference/ensemble_solution_gpu.md)
  is retained as a deprecated alias that warns and forwards, and will be
  removed in a future release.
- New exported helper
  [`pestpp_available()`](https://max578.github.io/PESTO/reference/pestpp_available.md)
  – a non-erroring probe for a PEST++ family executable
  (e.g. `pestpp-ies`, `pestpp-glm`). It is the documented way for
  examples, vignettes, and conditional tests to skip gracefully when no
  external binary is installed; every PESTO algorithm runs natively in R
  without one.
- The comparison vignette is reframed as *PEST and PEST++ Comparison and
  Simulation Study* and gains a head-to-head **cross-tool benchmark**
  (PESTO vs classic PEST 18.25 vs `pestpp-ies` 5.2.16) on a well-posed
  linear problem and a non-linear ODE. The figures are frozen real
  outputs of a fixed-seed reproducibility harness, shipped in
  `inst/extdata/pestpp_cache/`: accuracy parity on the linear problem,
  the ensemble methods’ advantage over linearised GLM on the non-linear
  one, a two-to-three-orders-of-magnitude wall-clock advantage, and the
  calibration caveat (raw ensemble intervals under-cover; apply
  inflation). A new *Lineage and scope* section grounds PESTO in PEST
  (Doherty 2015) and PEST++ (White et al. 2020) and states the algorithm
  boundary that makes the comparison fair.
- Corrected the benchmarked `pestpp-ies` version string in the vignette
  to 5.2.16 (was mislabelled 5.2.25).
- **APSIM positioned as PESTO’s flagship simulator partner** across the
  public surface, without narrowing the model-independent identity. The
  README, the *Getting started* vignette, and the package DESCRIPTION
  now name APSIM (coupled in-process via
  [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  and the `apsimx` package) as the primary worked partnership, alongside
  hydrological models, other crop models, and ODE systems. The APSIM
  vignette is retitled *Calibrating APSIM with PESTO* and was previously
  missing from the README documentation list (now fixed). No code or API
  change.
- **New publication-grade APSIM case study.** A synthetic-truth recovery
  experiment (OSSE) calibrating real APSIM Wheat parameters with the
  IES: it recovers two strong, mechanistically distinct parameters (soil
  runoff and fertiliser nitrogen) with the truth inside the 90% credible
  interval, correctly flags a third (sowing depth) as weakly
  identifiable, demonstrates out-of-sample under-coverage and its
  correction by covariance inflation, and emits a
  `pesto_ensemble_manifest`. Shipped as the *Calibrating APSIM Wheat
  with PESTO* vignette (frozen real outputs, so it builds without APSIM)
  plus a path-free reproducible driver at
  `system.file("case_studies/apsim_wheat_calibration.R", package = "PESTO")`.
  A second part calibrates physiological parameters (RUE, phenology) to
  the **real observed biomass** in
  [`apsimx::obsWheat`](https://rdrr.io/pkg/apsimx/man/obsWheat.html):
  PESTO’s posterior brackets the independent
  [`apsimx::optim_apsimx()`](https://rdrr.io/pkg/apsimx/man/optim_apsimx.html)
  optimum for both parameters and predicts held-out dates, adding the
  uncertainty the point optimiser lacks (driver
  `apsim_wheat_realdata.R`).

## PESTO 0.8.0

### New features

- [`train_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/train_gp_surrogate_tuned.md)
  and its companion
  [`predict_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate_tuned.md)
  fit a GP surrogate at **maximum-likelihood length scales** rather than
  the median heuristic
  [`train_gp_surrogate()`](https://max578.github.io/PESTO/reference/train_gp_surrogate.md)
  uses by default. The fit is **anisotropic** by default – one length
  scale per input dimension, estimated by maximising the GP’s own log
  marginal likelihood through per-axis coordinate pre-scaling (no change
  to the C++ kernel), with the response centred first. On a strongly
  anisotropic response this is a large accuracy gain over a single
  length scale: on the Branin function it cuts held-out error about
  23-fold versus the median-heuristic default (to roughly half a percent
  of the function range) and about 4.6-fold versus a single MLE length
  scale, bringing the surrogate within a small factor of a dedicated
  anisotropic GP (`DiceKriging`). This resolves a known finding that the
  *default* surrogate was several-fold worse than an MLE GP on an
  anisotropic function. The GP remains zero-mean, so it stays a small
  factor above a fully-optimised trend-bearing GP – a documented
  limitation, not a length-scale defect.

## PESTO 0.7.0

The release that readies PESTO as the authoritative ensemble-manifest
emitter for coordinated multi-tool runs and broadens the forward-model
surface beyond the `apsimx` simulator to native
ordinary-differential-equation models.

### New features

- **ODE / compartmental forward-model templates.** New exported builders
  turn a system of ordinary differential equations into a typed
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  that plugs straight into
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md),
  a
  [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  stack, and the manifest emitter – the ODE analogue of the
  [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  adapter.
  [`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
  is the generic builder (supply a `function(t, y, theta)` right-hand
  side, the initial state, and the time grid);
  [`crop_growth_forward_model()`](https://max578.github.io/PESTO/reference/crop_growth_forward_model.md)
  is the logistic dry-matter-accumulation crop template (Goudriaan &
  Monteith 1990); and
  [`seir_forward_model()`](https://max578.github.io/PESTO/reference/seir_forward_model.md)
  is the closed-population SEIR epidemic template (Anderson & May 1991).
  Integration is a self-contained fixed-step RK4 by default (no new hard
  dependency); `solver = "desolve"` delegates to the optional `deSolve`
  package for stiff systems. Each template is exercised by a
  simulate-forward-then-invert test that recovers the generating
  parameters.

- **Manifest schema `1.1.0`: grounded semantic descriptor
  (`obs_schema`).** `pesto_ensemble_manifest` gains an optional
  `obs_schema` slot stating the physical quantity and unit of each
  output and parameter column (plus optional per-column grounding
  provenance: `verified_on`, `oracle_kind`, `evidence_path`). Build one
  with the new exported
  [`pesto_obs_schema()`](https://max578.github.io/PESTO/reference/pesto_obs_schema.md)
  and pass it through `as_manifest(fit, obs_schema = ...)`. The
  descriptor turns column meaning from out-of-band roxygen convention
  into a machine-checkable field, so a downstream consumer can verify
  two manifests are commensurable by name rather than positionally. The
  class validator rejects a descriptor naming a column that does not
  exist in the data. `obs_schema` is provenance metadata and is not
  folded into `data_hash` (correspondence is grounded by an independent
  consumer, not by self-hash). Additive and backward-compatible: a
  `1.0.0` manifest reads back unchanged with `obs_schema = NULL`.

- **[`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  stamps the in-use APSIM version.** The returned forward-model closure
  now carries an `"apsim_version"` attribute (from
  [`apsimx::apsim_version()`](https://rdrr.io/pkg/apsimx/man/apsim_version.html),
  `NA_character_` when undeterminable), so a calibrated run can be
  grounded to the exact simulator that produced it via
  `as_manifest(fit, apsim_version = attr(fm, "apsim_version"))`.

### Minor improvements and fixes

- The *In-Process IES via R Callback* vignette gains an
  over-determination guard section: it shows how conditioning on a
  likelihood tighter than the data deserve (passing the standard error
  of a replicate mean, \sigma/\sqrt{m}, instead of the field-realistic
  replicate spread \sigma) collapses the ensemble and produces a
  confidently-wrong posterior, and how
  [`ensemble_spread_ess()`](https://max578.github.io/PESTO/reference/ensemble_spread_ess.md)
  and credible-interval coverage diagnose it.

- [`pesto_obs_schema()`](https://max578.github.io/PESTO/reference/pesto_obs_schema.md)
  gains a runnable example.

## PESTO 0.6.0

### Covariance inflation and localisation against ensemble collapse

Adds two opt-in countermeasures to the finite-ensemble pathologies that
make an iterative ensemble smoother over-confident: covariance inflation
(against under-dispersion / ensemble collapse) and covariance
localisation (against spurious finite-sample parameter-observation
correlations). Both default to off; a `NULL` specification leaves
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
and
[`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md)
byte-identical to the previous release.

#### New exports

- [`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)
  – inflation specification with four methods: `"rtps"` (relaxation to
  prior spread, Whitaker & Hamill 2012; the per-parameter,
  spectrally-aware workhorse), `"adaptive"` (global inflation targeting
  a spread-retention floor), `"multiplicative"` (fixed factor), and
  `"none"`.
- [`pesto_localisation()`](https://max578.github.io/PESTO/reference/pesto_localisation.md)
  – localisation specification: `"correlation"` (automatic,
  coordinate-free, Luo & Bhakta 2020 – the recommended default for
  parameter problems with no spatial metric) or `"distance"` (classical
  Gaspari-Cohn taper of a parameter-to-observation distance matrix).
- [`ensemble_spread_ess()`](https://max578.github.io/PESTO/reference/ensemble_spread_ess.md)
  – the collapse *diagnostic*: the spectral participation ratio of the
  parameter anomaly covariance, i.e. the effective number of
  variance-carrying directions. Recorded on every iteration regardless
  of method.
- [`correlation_localisation()`](https://max578.github.io/PESTO/reference/correlation_localisation.md),
  [`gaspari_cohn()`](https://max578.github.io/PESTO/reference/gaspari_cohn.md),
  [`ensemble_solution_localised()`](https://max578.github.io/PESTO/reference/ensemble_solution_localised.md)
  – the C++ kernels backing the above.
  [`ensemble_solution_localised()`](https://max578.github.io/PESTO/reference/ensemble_solution_localised.md)
  is the explicit-gain GLM update that hosts the Schur-product
  localisation the SVD kernel cannot; with no taper it reproduces
  [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  (approximate form) to truncation tolerance.

#### Behaviour

- [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  and
  [`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md)
  gain `inflation` and `localisation` arguments (both `NULL` by default)
  and now record the spread-ESS and (when active) inflation /
  localisation diagnostics in their per-step metadata, which flow into
  the ensemble manifest.
- Localisation uses the approximate (upgrade_1) update; combining a
  non-`NULL` `localisation` with `use_approx = FALSE` warns and drops
  the null-space correction.

Note on terminology: the spectral spread-ESS is scale-invariant, so it
is used as the collapse *diagnostic*, while the `"adaptive"` inflation
targets a variance-*magnitude* retention floor; `"rtps"` is the method
that reshapes the spectrum. See the *Countering Ensemble Collapse:
Inflation and Localisation* vignette.

## PESTO 0.5.0

### Forward-model contract + first-class multi-fidelity bridge

Promotes the two-adapter forward-model contract from an implicit
convention to a typed, enforceable object, and makes the multi-fidelity
`(cheap, expensive)` bridge first-class (APSIM-bridge invariants 1 and
3). No breaking changes to existing calls:
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
still accepts a bare `function(theta) -> obs`.

#### New exports

- [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  — an S7 class wrapping a forward callable with its output
  dimensionality, expected parameter names, failure policy
  (`on_failure`, `max_fail_frac`), evaluation strategy (serial /
  `"multicore"` / custom `map_fn`), and a `fidelity` tag. This is the
  single contract both the native-callback and `.pst`-file adapter modes
  honour.
- [`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md)
  — generic that runs a forward model (or a multi-fidelity model at a
  chosen `level`) and returns a shape-guaranteed `nreal x nobs` matrix
  with `"n_failures"` / `"fail_idx"` attributes.
- [`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md)
  — coerces a bare function (or passes through an existing object) into
  the contract; used internally so bare functions keep working
  unchanged.
- [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  — an ordered stack of fidelity levels (cheapest first) plus relative
  `costs`; the first-class form of the bridge’s fidelity vector.
- [`mf_control_variate()`](https://max578.github.io/PESTO/reference/mf_control_variate.md)
  — the affine (Kennedy-O’Hagan AR(1)) control-variate primitive that
  debiases a cheap level against a sparse expensive sample; the plug-in
  point for surrogate cascades.

#### Sequential (filter-mode) IES

- New
  [`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md)
  — a filtering counterpart to the batch smoother
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).
  It assimilates time-ordered observation `windows` one after another
  against a static parameter ensemble, the posterior of each window
  becoming the prior of the next, so a tightening parameter posterior is
  available after every window (the in-season assimilation case). It
  reuses the forward-model contract (parallel- and multi-fidelity-ready
  via a per-window `fidelity_schedule`) and the C++
  [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  kernel; `window_noptmax > 1` gives an iterated filter per window. The
  result records a per-window history including the per-parameter
  ensemble standard deviation (the tightening trace).
- Filter results (`pesto_ies_filter_result`) flow into the manifest
  contract:
  [`as_manifest()`](https://max578.github.io/PESTO/reference/as_manifest.md)
  tags them `method = "ies_filter"` (added to the
  `pesto_ensemble_manifest` validator) and carries their fidelity
  provenance, so a filtered ensemble is a first-class scenario for any
  downstream consumer.

#### Behaviour

- [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  gains `fidelity_schedule` (consulted only for a
  `pesto_multifidelity_model`): the fidelity level evaluated at each
  iteration, supporting cheap-early / expensive-late ramping. The final
  ensemble refresh always uses the highest fidelity.
- Fidelity provenance now closes the manifest (C2) lineage: a
  multi-fidelity
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  run records its realised schedule in the result
  (`$fidelity = list(type, schedule, final_level, n_levels, costs)`),
  [`as_manifest()`](https://max578.github.io/PESTO/reference/as_manifest.md)
  inherits it into the `pesto_ensemble_manifest` `fidelity` slot unless
  overridden, and
  [`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md)
  /
  [`read_manifest()`](https://max578.github.io/PESTO/reference/read_manifest.md)
  round-trip the structured record faithfully (it is outside the
  integrity hash, so it does not affect
  [`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)).
  Single-fidelity runs record `NULL`, so their manifests are unchanged.
  The manifest `fidelity` slot is now documented as a structured
  provenance list (legacy named-numeric tags are still accepted on
  read).
- Parallel, fault-tolerant ensemble evaluation: a `pesto_forward_model`
  with `parallel = "multicore"` dispatches realisations across forked
  workers via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  with L’Ecuyer streams (reproducible under `RNGkind("L'Ecuyer-CMRG")`);
  serial bulk evaluation is unchanged and remains the default.
- [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  now writes each realisation to a unique per-run file, making the
  closure safe to drive in parallel (wrap it in a
  `pesto_forward_model(parallel = "multicore")`).
- The internal bulk-then-per-row evaluation helper
  (`.eval_forward_safe`) was retired in favour of the shared engine
  behind
  [`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md);
  the on-error abort message changed from `` `forward_model` failed ``
  to `forward model failed`.

#### Dependencies

- `parallel` (a base R package) added to `Imports`.

## PESTO 0.4.1

### AAGI guidance uplift

A code-aesthetics and review-readability patch on top of `0.4.0`. No
runtime behaviour changes; no exported-API changes; no shipped-data
changes. The aim is to lift every source surface to the bar set out in
`r_style.md` ahead of any AAGI-AUS push.

#### Validation delegation (r_style.md invariant 9)

- `R/internal_validation.R` introduces the shared primitive validators
  (`.assert_positive_scalar()`, `.assert_nonneg_scalar()`,
  `.assert_character_scalar()`, `.assert_logical_scalar()`,
  `.assert_path_exists()`, `.assert_matrix()`,
  `.assert_numeric_vector()`, `.assert_function()`,
  `.assert_data_frame()`, `.assert_choice()`, `.assert_same_ncol()`,
  `.assert_same_nrow()`, `.assert_required_cols()`). All
  `@noRd `[`@keywords`](https://github.com/keywords)` internal`; every
  helper signals failure via `stop(call. = FALSE, ...)` with a
  backticked argument name.
- Public functions in `apsim_callback.R`, `pesto_reference_ies.R`,
  `pesto_run.R` (`pesto_ies_callback`), `pst_io.R`, `scenario.R`,
  `surrogate.R`, `manifest.R`, `ensemble_io.R`, `plot.R`, and
  `check_surrogate_regime.R` now open with `.check_*` / `.assert_*`
  calls instead of inline `if (!is.x) stop(...)` walls.
- Error messages backtick the offending argument name throughout (Sparks
  convention).

#### Section banners (r_style.md invariant 4)

- Long function bodies (`pesto_ies_callback`, `pesto_ies`, `pesto_glm`,
  `pesto_sweep`, `pesto_sensitivity`, `read_pst`, `write_pst`,
  `apsim_callback`, `pesto_reference_ies`, `.find_pestpp_exe`) now carry
  Sparks-style dash-banner section comments that paragraph the work
  (validate inputs / resolve paths / iterate / parse outputs / assemble
  result).

#### Import concentration (r_style.md invariant 10)

- `@importFrom ggplot2` annotations consolidated into
  `R/pesto-package.R`; the per-function `@importFrom` annotation on
  [`plot_phi()`](https://max578.github.io/PESTO/reference/plot_phi.md)
  has been removed in favour of inline `ggplot2::` qualification at the
  call sites.

#### Vignette prose

- Prose semicolons in `vignettes/apsim-callback.Rmd` and
  `vignettes/ensemble-manifest.Rmd` converted to `. Capital` joins per
  `manuscript_style.md` invariant 5.

#### README

- Added `Dependencies`, `Contributing`, and `Acknowledgements` sections
  per the AAGI repository-guidelines README contract. Citation block
  bumped to `R package version 0.4.1`.

#### Canonical metadata URLs

- `DESCRIPTION URL` and `BugReports`, `CITATION.cff`, `codemeta.json`,
  `inst/CITATION`, `_pkgdown.yml`, and the README citation + issues
  links now point to `https://github.com/AAGI-AUS/PESTO` (canon
  checklist item 5). README install instructions and the personal
  r-universe URL are retained at `max578/PESTO` and
  `https://max578.r-universe.dev` as interim distribution infrastructure
  until the AAGI-AUS push lands.

## PESTO 0.4.0

### AAGI recipes uplift and canon channel migration

This release contains no R, C++, or shipped-data changes. It is a
governance, metadata, and project-hygiene release that lands the AAGI
canon recipes on the `max578/PESTO` channel.

#### Canon channel migration

- Primary canon channel migrated from `AAGI-AUS/PESTO` to
  `max578/PESTO`. `DESCRIPTION`, `CITATION.cff`, `codemeta.json`,
  `inst/CITATION`, `_pkgdown.yml`, `README.md`, `CONTRIBUTING.md`,
  `API_STABILITY.md`, and the pkgdown GitHub Actions workflow header now
  point to `https://github.com/max578/PESTO` and
  `https://max578.github.io/PESTO`. The `aagi` git remote is retained as
  a frozen read-only mirror; no push to `AAGI-AUS` without explicit
  per-instance maintainer approval.
- The package opts out of the `AAGI-AUS` publication canon (it is
  published through the `max578` channel); project-local configuration
  is excluded from R-package builds via `.Rbuildignore`.
- `man/PESTO-package.Rd` regenerated to inherit the new URLs from
  `DESCRIPTION` via `devtools::document()`.

#### Metadata version sync

- `CITATION.cff` (`version: 0.1.0`), `codemeta.json`
  (`version: "0.1.0"`), `inst/CITATION` (`R package version 0.3.3`), and
  the README citation block were not in lock-step with `DESCRIPTION`.
  All four are now on `0.4.0` with `date-released: "2026-05-28"` and
  `dateModified: "2026-05-28"`.

#### Sole copyright holder

- `codemeta.json` `copyrightHolder` corrected to `Person "Max Moldovan"`
  with ORCID `0000-0001-9680-8474` and Adelaide University affiliation,
  matching `Authors@R` and `LICENSE.md`.

#### AAGI canon recipes added

- `CODE_OF_CONDUCT.md` (Contributor Covenant v2.1, pointer form).
- `SECURITY.md` (vulnerability reporting policy; maintainer email,
  five-working-day acknowledgement, scope statement).
- `air.toml` (Air formatter configuration: 80-char line width, two-space
  indent, auto line endings).
- `.lintr` (lintr defaults aligned with `r_style.md` direction: 80-char
  line, `snake_case` / `dotted.case` / symbols object names, two-space
  indent; `src`, `tools`, `inst/extdata`, `vignettes` excluded).
- `.Rbuildignore` already excluded all four paths; no tarball impact.

## PESTO 0.3.3

### FLIBS Makevars portability fix (closes critical-review P2 [\#7](https://github.com/max578/PESTO/issues/7))

- `src/Makevars` `PKG_LIBS` now follows `$(BLAS_LIBS)` with `$(FLIBS)`
  per Writing R Extensions §1.2.1.5. Resolves the pre-existing
  structural WARNING (“apparently using \$(BLAS_LIBS) without following
  \$(FLIBS) in ‘src/Makevars’”) that had survived several check passes
  as a “documented local-env artifact”. `Makevars.win` was already
  correct; no change needed there.
- The remaining baseline structural WARNING (Apple-clang 21.0.0 vs.
  `R_ext/Boolean.h:62` `-Wfixed-enum-extension` pragma) is in R’s own
  header on this toolchain version and is not present on CRAN’s build
  farm; it persists harmlessly.

## PESTO 0.3.2

### CSV-only manifest mode renamed (post critical-review, 2026-05-16)

- `write_manifest(format = "csv")` is renamed to
  `write_manifest(format = "csv_unverified")` to flag the weaker
  integrity contract at every call-site. The mode itself is unchanged
  (CSV-only sidecars, hash recorded but not disk-verifiable).
- New top-level YAML field `integrity: verifiable | not_verifiable`
  derived from `format`. Verifiable: `rds`, `both`. Not verifiable:
  `csv_unverified`. Lets non-R downstream tools (e.g. Python pipelines)
  branch on the integrity contract without parsing the PESTO-specific
  format vocabulary.
- The legacy `format = "csv"` spelling is still accepted at the API
  boundary with a deprecation warning; the persisted form always uses
  `csv_unverified`.
  [`read_manifest()`](https://max578.github.io/PESTO/reference/read_manifest.md)
  normalises old YAMLs on the read side, so 0.3.1 manifests round-trip
  cleanly under 0.3.2.
- Validator vocabulary updated: `{rds, both, csv_unverified}`. The old
  `"csv"` token is rejected at slot-set time (only the renamed argument
  accepts it, with a warning).
- New tests cover the renamed value, the deprecation-warning path, and
  the `integrity:` YAML field for verifiable modes.
- Vignette `ensemble-manifest.Rmd` reframed: explicit “Inspection CSVs
  (verifiable, via format = ‘both’)” vs. “Unverified CSV export (via
  format = ‘csv_unverified’)” sections; the latter is presented as “for
  export, not for storage you intend to re-load and trust”.

## PESTO 0.3.1

### Manifest sidecar `format=` option (roadmap §A5 polish)

- [`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md)
  gains a `format = c("rds", "both", "csv")` argument. `"rds"` (default)
  preserves the current bit-exact binary behaviour. `"both"` writes RDS
  sidecars plus parallel CSV inspection files (`*_inspection.csv`); the
  SHA-256 hash stays bound to the RDS. `"csv"` writes CSV-only sidecars
  for inspection / interchange workflows where bit-exact integrity is
  not required.
- New S7 slot `format` on `pesto_ensemble_manifest` records the on-disk
  serialisation mode (default `"rds"`; preserved through read/write
  round-trips). Validator enforces the three-value vocabulary.
- [`read_manifest()`](https://max578.github.io/PESTO/reference/read_manifest.md)
  dispatches on file extension in the YAML’s `artefacts:` block — reads
  RDS via [`readRDS()`](https://rdrr.io/r/base/readRDS.html), CSV via
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html).
- [`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
  gains a `message` field on its return list and returns `ok = NA` (with
  explanation) for `format = "csv"` manifests whose IEEE 754 doubles
  have round-tripped through a write formatter. Existing
  `format = "rds"` callers see no behaviour change; the new field is
  `NULL` in that case.
- YAML schema picks up a top-level `format:` key and an optional
  `inspection_csv:` block when `format = "both"`. Backwards-compatible:
  YAMLs written by PESTO 0.3.0 (no `format:` key) read back with
  `format = "rds"` per the default. No schema-version bump required.
- New tests in `test-manifest.R` cover all three formats plus the
  unknown-format rejection path.

## PESTO 0.3.0

### Ensemble Manifest as a Portable, Versioned Run Record

- New S7 class `pesto_ensemble_manifest` — versioned, hashed,
  provenance-tracked container for ensemble-run output. Slots cover
  `params`, `outputs`, `weights`, `obs_target`, `seed`, `data_hash`
  (SHA-256), `fidelity`, `apsim_version`, `pesto_version`, `timestamp`,
  plus method context (`method`, `noptmax`, `lambda_schedule`,
  `failure_rate`). This is the documented, versioned format that any
  downstream tool can read.
- New
  [`as_manifest()`](https://max578.github.io/PESTO/reference/as_manifest.md)
  — S7 generic with a method for `pesto_ies_callback_result`.
  Non-destructive: wraps without mutating the source result.
- New
  [`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md)
  /
  [`read_manifest()`](https://max578.github.io/PESTO/reference/read_manifest.md)
  — YAML+RDS serialisation. The YAML carries metadata + relative paths
  to three sidecar RDS files (`*_params.rds`, `*_outputs.rds`,
  `*_assim.rds`); RDS is used in preference to CSV so IEEE 754 doubles
  round-trip bit-exactly (the SHA-256 integrity check would otherwise
  trip on CSV formatter precision loss).
- New
  [`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
  — recomputes the SHA-256 over
  `(params, outputs, weights, obs_target, seed)` and compares to the
  stored value, returning a diagnostic list. Detects post-write
  tampering with the sidecar CSVs.
- New vignette: `ensemble-manifest` — end-to-end demo of construct →
  write → read → verify, plus tamper-detection.

### In-Process IES via R Callback (UQ ag-stack roadmap §D4)

- New
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  — drives an Iterative Ensemble Smoother entirely in R using a
  user-supplied forward-model callable, bypassing the `.pst`-file
  write/read cycle of
  [`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md).
  Each iteration calls the existing C++ kernel
  [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  (Chen & Oliver, 2013). Tolerates per-realisation failures via
  `on_failure = c("na", "stop")` and reports a `failure_rate` in the
  result object. Phase-1 behaviour uses a single lambda per iteration
  (or user-supplied schedule); a full `pestpp-ies`-style lambda
  line-search is a planned Phase-2 enhancement.
- New
  [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  — adapter that wraps the `apsimx` package (now in `Suggests`) into a
  forward-model closure suitable for
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).
  Per-realisation template copy, parameter edit via
  [`apsimx::edit_apsimx()`](https://rdrr.io/pkg/apsimx/man/edit_apsimx.html)
  / [`edit_apsim()`](https://rdrr.io/pkg/apsimx/man/edit_apsim.html),
  run, and extraction. Failures (edit / run / extractor) surface as `NA`
  rows for the IES driver to handle.
- New vignette: `apsim-callback` — synthetic linear-Gaussian recovery
  demo plus disabled apsimx example.
- New benchmark script: `inst/benchmarks/d4_callback_vs_pst.R`.
- [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  records `obs_target`, `obs_sd`, and `weights` on its result list so
  the manifest emitter has full IES context to capture.

### Notes

- Real-APSIM benchmarking (target ≥10× speed-up vs `.pst` path) is
  deferred to the §D1 scenario library landing. The current benchmark
  script measures the callback path on a synthetic surrogate forward
  model only.
- `apsimx` is in `Suggests` not `Imports`;
  [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  checks for it at call time with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html).
- New hard imports: `S7 (>= 0.2.0)`, `yaml (>= 2.3.0)`,
  `digest (>= 0.6.0)`.

## PESTO 0.2.0 (2026-04-25)

### Bug fixes

- **I1 —
  [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  sign-convention bug.** The C++ kernel requires
  `obs_resid = sim - obs`; the docstring previously stated the inverse.
  Two genuine in-package call sites were silently inverting upgrades:
  `src/surrogate_ies.cpp:347` and
  `vignettes/surrogate-ies.Rmd:148, 181`. Both fixed; surrogate-IES now
  applies upgrades in the correct direction. Regression test
  `tests/testthat/test-ensemble-solution-sign.R` asserts strict monotone
  phi descent under the correct convention AND geometric divergence
  under the inverted one.

### New exported functions

- [`pesto_reference_ies()`](https://max578.github.io/PESTO/reference/pesto_reference_ies.md)
  — pure-R, textbook implementation of the Chen & Oliver (2013) eq. 12
  IES update. Independent of the C++ kernel; used as the canonical
  comparison target by the comparison vignette so it ships and runs
  without the upstream `pestpp-ies` binary. Cross-validated against the
  C++ kernel at machine precision (max element-wise delta = 5.8e-15).
- [`check_surrogate_regime()`](https://max578.github.io/PESTO/reference/check_surrogate_regime.md)
  — soft guardrail that warns when the surrogate-IES regime is
  unfavourable (`n_train < threshold * n_params`). Stand-alone helper,
  not auto-invoked by
  [`pesto_surrogate_ies()`](https://max578.github.io/PESTO/reference/pesto_surrogate_ies.md)
  (v0.3 wiring candidate).

### API enhancements

- [`plot_identifiability()`](https://max578.github.io/PESTO/reference/plot_identifiability.md)
  gains a `jacobian = NULL` matrix-input path. Backward-compatible:
  `jco_file = NULL` retained; the two are mutually exclusive.

### Self-contained `pestpp-ies` comparison

- `vignettes/pestpp-comparison-and-simulation.Rmd` now compares PESTO
  native IES against
  [`pesto_reference_ies()`](https://max578.github.io/PESTO/reference/pesto_reference_ies.md)
  by default — no upstream binary required. The pure-R reference cache
  ships at `inst/extdata/pestpp_cache/scenario_a_reference.rds`
  (SHA-256-pinned to the prior ensemble). When the developer-side cache
  `tools/pestpp_benchmark/scenario_a_pestpp_ies.rds` is present and
  `PESTO_PESTPP_BIN` resolves, the vignette extends the agreement plot
  with the live binary’s posterior.
- Hardcoded local absolute path replaced with
  `Sys.getenv("PESTO_PESTPP_BIN")` + `Sys.which("pestpp-ies")` fallback.
- New `tools/pestpp_benchmark/run_benchmark.R` regenerates both caches
  deterministically. Documented in `CONTRIBUTING.md`.

### Documentation and CRAN-readiness

- Every export now has a runnable `@examples` block (30 of 30 documented
  exports/methods). The four external-binary runners (`pesto_ies`,
  `pesto_glm`, `pesto_sweep`, `pesto_sensitivity`) and
  `pesto_surrogate_ies` use guarded `\donttest{}` (no `\dontrun{}`).
- Vignettes acquire a “Regime of applicability” subsection
  (`surrogate-ies.Rmd`) and a “surrogate savings in this regime” note
  (`pestpp-comparison-and-simulation.Rmd`) covering the
  curse-of-dimensionality finding from investigation I3.
- Kernel docstring
  [`?ensemble_solution`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  now states the `sim - obs` convention with a full GLM-derivation
  rationale.
- New `cran-comments.md` with per-NOTE justification.
- New `CITATION.cff` (CFF 1.2.0 + ORCID + preferred citation).
- New `codemeta.json` (CodeMeta 2.0).
- New `CONTRIBUTING.md` documenting the developer benchmark workflow.
- New `inst/WORDLIST` with ~120 domain terms; `Language: en-AU` added to
  DESCRIPTION.

### Build and CRAN-portability

- `src/Makevars`: `PKG_LIBS` gains `$(FLIBS)` (CRAN portability
  requirement for `$(BLAS_LIBS)`).
- `LICENSE` renamed to `LICENSE.md` and `.Rbuildignore`-d (CRAN
  convention for GPL-licensed packages).

### Open investigations (status)

- I1 — **CLOSED** (this release).
- I2 — partially closed: PESTO C++ kernel verified equivalent to Chen &
  Oliver (2013) eq. 12 at machine precision. Residual
  PESTO-vs-pestpp-ies posterior gap is structural (Marquardt
  sub-cycling + per-realisation noise perturbation in upstream
  pestpp-ies) and is a v0.3 implementation milestone.
- I3 — **DOCUMENTED** (this release): regime-of-applicability text added
  to both vignettes;
  [`check_surrogate_regime()`](https://max578.github.io/PESTO/reference/check_surrogate_regime.md)
  helper exported.

## PESTO 0.1.1

### Sole-Authorship Consolidation

- Authorship consolidated to Max Moldovan as sole `aut`, `cre`, and
  `cph` (administrative consolidation; no licence change).
- Licence unchanged: GPL-3 or any later version.
- `LICENSE` file rewritten with corrected canonical wording and
  copyright attribution.
- Source-file headers in `src/` updated to reflect sole authorship.

## PESTO 0.1.0

### Initial Release

#### Core Features

- [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  — High-performance C++ implementation of the IES ensemble update
  equation (Chen & Oliver, 2013) via RcppEigen.
- [`ensemble_solution_mda()`](https://max578.github.io/PESTO/reference/ensemble_solution_mda.md)
  — Multiple Data Assimilation (Evensen, 2018) update kernel.
- [`compute_phi()`](https://max578.github.io/PESTO/reference/compute_phi.md)
  — Fast weighted sum-of-squares objective function.

#### Hardware Acceleration

- [`adaptive_svd()`](https://max578.github.io/PESTO/reference/adaptive_svd.md)
  — Automatic SVD backend selection (LAPACK, Eigen BDCSVD, or randomised
  SVD) based on matrix size and target rank.
- [`rsvd()`](https://max578.github.io/PESTO/reference/rsvd.md) —
  Randomised SVD (Halko-Martinsson-Tropp, 2011) for asymptotically
  faster rank-k approximations.
- [`accelerate_svd()`](https://max578.github.io/PESTO/reference/accelerate_svd.md)
  — Direct LAPACK SVD leveraging platform-optimised BLAS (Apple
  Accelerate/AMX on macOS, MKL or OpenBLAS on Linux).
- [`ensemble_solution_gpu()`](https://max578.github.io/PESTO/reference/ensemble_solution_gpu.md)
  — ensemble solution with adaptive SVD backend selection and
  performance diagnostics (renamed to
  [`ensemble_solution_adaptive()`](https://max578.github.io/PESTO/reference/ensemble_solution_adaptive.md)
  in the development version).

#### Surrogate-Accelerated IES (Novel)

- [`train_gp_surrogate()`](https://max578.github.io/PESTO/reference/train_gp_surrogate.md)
  — Gaussian Process surrogate model training with automatic
  hyperparameter selection (median heuristic).
- [`predict_gp_surrogate()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate.md)
  — GP prediction with uncertainty quantification.
- [`surrogate_ensemble_update()`](https://max578.github.io/PESTO/reference/surrogate_ensemble_update.md)
  — Surrogate-accelerated IES update with adaptive model/surrogate
  switching and control-variate bias correction.
- [`adaptive_ensemble_size()`](https://max578.github.io/PESTO/reference/adaptive_ensemble_size.md)
  — Convergence-aware dynamic ensemble sizing based on ESS and
  coefficient of variation diagnostics.

#### PEST++ Integration

- [`read_pst()`](https://max578.github.io/PESTO/reference/read_pst.md) /
  [`write_pst()`](https://max578.github.io/PESTO/reference/write_pst.md)
  — PEST control file I/O.
- [`read_ensemble()`](https://max578.github.io/PESTO/reference/read_ensemble.md)
  /
  [`write_ensemble()`](https://max578.github.io/PESTO/reference/write_ensemble.md)
  — Ensemble file I/O (CSV + binary).
- [`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md),
  [`pesto_glm()`](https://max578.github.io/PESTO/reference/pesto_glm.md),
  [`pesto_sweep()`](https://max578.github.io/PESTO/reference/pesto_sweep.md),
  [`pesto_sensitivity()`](https://max578.github.io/PESTO/reference/pesto_sensitivity.md)
  — High-level wrappers for PEST++ executables.
- [`create_pest_scenario()`](https://max578.github.io/PESTO/reference/create_pest_scenario.md)
  — Programmatic scenario builder.

#### Visualisation

- [`plot_phi()`](https://max578.github.io/PESTO/reference/plot_phi.md) —
  Objective function convergence plotting.
- [`plot_ensemble()`](https://max578.github.io/PESTO/reference/plot_ensemble.md)
  — Prior/posterior parameter distribution comparison.
- [`plot_identifiability()`](https://max578.github.io/PESTO/reference/plot_identifiability.md)
  — SVD-based parameter identifiability analysis.
- [`plot_surrogate_diagnostics()`](https://max578.github.io/PESTO/reference/plot_surrogate_diagnostics.md)
  — Surrogate IES performance visualisation.
