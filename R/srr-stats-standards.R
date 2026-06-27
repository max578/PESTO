#' srr_stats
#'
#' Compliance documentation for the rOpenSci statistical-software standards.
#' `@srrstats` tags lead with how PESTO complies; non-applicable standards are
#' in the NA_standards block. Some standards are documented next to the code or
#' test that satisfies them (see the breadcrumbs).
#'
#' @srrstatsVerbose TRUE
#' @srrstats {G1.0} Primary references are cited in roxygen (Chen & Oliver 2013; Evensen 2018; White 2018; Goudriaan & Monteith 1990; Anderson & May 1991).
#' @srrstats {G1.1} PESTO is the first R-native implementation of PEST++-class iterative ensemble smoothers (previously available only in C++/PEST++), adding original extensions (typed forward-model contract, in-process callback, multi-fidelity, surrogate acceleration); stated in DESCRIPTION and README.
#' @srrstats {G1.2} A Life Cycle Statement is provided in API_STABILITY.md (additive-then-lifecycle policy).
#' @srrstats {G1.3} Statistical terms (IES, MDA, objective function phi, prior/posterior ensemble, inflation, localisation, spread-ESS) are defined where introduced in the function documentation and vignettes.
#' @srrstats {G1.4} All exported functions are documented with roxygen2 (see man/).
#' @srrstats {G1.4a} Internal (non-exported) helpers carry @noRd and @keywords internal (e.g. R/internal_validation.R).
#' @srrstats {G1.5} The performance claims (vs PEST 18.25 / pestpp-ies 5.2.16) are reproduced in the pestpp-comparison-and-simulation vignette and an in-suite stored-golden test (test-correctness-analytic.R).
#' @srrstats {G2.0} Input lengths are asserted via the .assert_*_scalar helpers (length(x) != 1L) in R/internal_validation.R.
#' @srrstats {G2.0a} Length expectations are documented in @param (e.g. obs_sd is a scalar or length-nobs vector; parcov is length-npar).
#' @srrstats {G2.1} Input types are asserted (numeric / character / logical / function / matrix) in R/internal_validation.R.
#' @srrstats {G2.1a} Data-type expectations are documented in @param (obs a named numeric vector; prior_ensemble a numeric matrix / data.table).
#' @srrstats {G2.2} Univariate parameters are guarded by scalar length assertions (R/internal_validation.R).
#' @srrstats {G2.3} Univariate character inputs are restricted via match.arg() / .assert_choice().
#' @srrstats {G2.3a} match.arg() is used for character-parameter dispatch (method, format, on_failure, taper).
#' @srrstats {G2.3b} Unmatched character values raise an explicit error via match.arg().
#' @srrstats {G2.4} Conversions are explicit where required (see G2.4a/G2.4b); factor/character conversions are not applicable (NA_standards block).
#' @srrstats {G2.4a} Explicit integer conversion via as.integer() (noptmax, fidelity levels).
#' @srrstats {G2.4b} Explicit numeric conversion via as.numeric() (obs, obs_sd, parcov).
#' @srrstats {G2.6} One-dimensional inputs are coerced uniformly: observations pass through as.numeric() regardless of their incoming 1-D class.
#' @srrstats {G2.7} The prior ensemble is accepted as a matrix, data.frame, or data.table (coerced to a numeric matrix in pre-processing).
#' @srrstats {G2.8} Inputs are converted to a single internal representation (a numeric matrix) during pre-processing before any sub-function or the C++ kernel sees them.
#' @srrstats {G2.10} Parameter-column extraction yields a numeric matrix consistently whether the prior is a matrix, data.frame, or data.table (the data.table path is tested in test-ies-callback.R).
#' @srrstats {G2.13} Missing data are checked in pre-processing: each iteration retains only all-finite realisations before the update (R/pesto_run.R), and input NA handling is governed by on_failure.
#' @srrstats {G2.14} Users control missing-data handling via on_failure: error or ignore (see G2.14a/G2.14b); imputation is not applicable (NA_standards block).
#' @srrstats {G2.14a} on_failure = "stop" errors on missing / failed realisations.
#' @srrstats {G2.14b} on_failure = "na" (default) carries failed realisations forward and records the failure_rate diagnostic.
#' @srrstats {G2.15} The driver never assumes non-missingness: it subsets to all-finite realisations before computing phi or passing data to the C++ update, so base reductions never see missing values.
#' @srrstats {G2.16} Undefined values are handled: non-finite forward outputs (NaN, Inf) are excluded each iteration via the all-finite realisation check, and Inf observation SDs (obs_sd = Inf -> zero weight) are supported (tested in test-bayesian-recovery.R).
#' @srrstats {G3.0} No naive floating-point equality is used (verified by source scan); tolerances / inequalities are used throughout.
#' @srrstats {G3.1} PESTO does not use stats::cov: covariances are formed from ensemble anomaly cross-products intrinsic to the IES algorithm, and the prior covariance diagonal is user-overridable via parcov.
#' @srrstats {G3.1a} The parcov override (the user-specifiable prior covariance diagonal) is documented at @param in ?pesto_ies_callback.
#' @srrstats {G4.0} Functions writing output files derive appropriate suffixes from the supplied path: write_manifest() builds .rds / .csv sidecars from the YAML basename (tools::file_path_sans_ext), and write_pst() writes .pst control files.
#' @srrstats {BS1.0} "hyperparameter" is used only for the GP-surrogate kernel parameters (length-scale, input-scale), defined in train_gp_surrogate_tuned() documentation.
#' @srrstats {BS1.1} Data entry (a prior parameter ensemble, a named observation vector, obs_sd) is shown in the pesto_ies_callback() examples, the README, and the getting-started vignette.
#' @srrstats {BS1.2} The prior is supplied as a prior parameter ensemble; how to specify it is documented at the prior_ensemble parameter with worked examples.
#' @srrstats {BS1.2a} Prior specification is described, with example code, in the README ("Ensemble Bayesian inference (IES)" section).
#' @srrstats {BS1.2b} The getting-started vignette shows prior vs posterior ensembles via plot_ensemble().
#' @srrstats {BS1.2c} pesto_ies_callback() examples construct a prior ensemble and run the smoother.
#' @srrstats {BS1.3} The computational-process parameters (noptmax, lambda, phi_tol, parcov, eigthresh, use_approx) are each documented in roxygen.
#' @srrstats {BS1.3a} Using a previous run's posterior ensemble as the starting (prior) ensemble of the next run (warm-start) is documented at @param prior_ensemble and tested (test-bayesian-interface.R).
#' @srrstats {BS1.3b} The alternative drivers (in-process callback, sequential filter, MDA, .pst path) are documented and cross-linked via @seealso.
#' @srrstats {BS2.1} Dimensional commensurability of inputs (nobs / npar / nreal) is enforced in the driver pre-processing.
#' @srrstats {BS2.1a} The dimension / length validation is tested (test-ies-callback.R, test-edge-conditions.R).
#' @srrstats {BS2.2} obs_sd (the observation-error specification) is validated as a distinct pre-processing step before the iteration loop.
#' @srrstats {BS2.3} obs_sd must be scalar or length-nobs; a non-conforming length raises an error (no excess values silently discarded).
#' @srrstats {BS2.4} obs_sd length is checked for commensurability with the observations (scalar recycled to nobs).
#' @srrstats {BS2.5} Variance-like parameters are constrained positive: obs_sd > 0 and parcov > 0 are enforced.
#' @srrstats {BS2.6} Computational parameters are range-checked (noptmax >= 1, positive eigthresh, positive phi_tol, in-range fidelity levels).
#' @srrstats {BS2.11} Starting points are supplied as the prior_ensemble matrix -- inherently plural, one row per realisation.
#' @srrstats {BS2.12} The verbose parameter controls progress output and defaults to TRUE.
#' @srrstats {BS2.15} on_failure = "na" catches forward-model failures and captures them in the return value (failure_rate + NA realisations) rather than aborting; tested in test-ies-callback.R.
#' @srrstats {BS3.0} Missing-value handling is documented via the on_failure parameter and the forward-model contract (failed realisations become NA rows, tolerated or aborted).
#' @srrstats {BS4.0} The assimilation algorithms are documented via citation (Chen & Oliver 2013; Evensen 2018) in roxygen references.
#' @srrstats {BS4.5} Non-convergence is surfaced via the returned phi trace, the per-iteration spread-ESS ratio, and the converged flag; failed forward evaluations are handled via on_failure; Marquardt lambda damping mitigates divergence.
#'
#' Standards documented next to their code/tests:
#' - G5 testing series: tests/testthat/ (correctness-analytic, recovery-scaling,
#'   edge-conditions, noise-susceptibility, helper-srr).
#' - BS2.7, BS2.8, BS2.13, BS5.1, BS5.3, BS5.5: test-bayesian-interface.R.
#' - BS4.1, BS7.2: test-correctness-analytic.R.
#' - BS4.2, BS7.0, BS7.1, BS7.4, BS7.4a: test-bayesian-recovery.R.
#' - BS7.3: test-recovery-scaling.R.
#' - BS6.0, BS6.1, BS6.3: R/ies_result_methods.R.
#' - BS1.4, BS4.3, BS4.4, BS4.6, BS4.7: test-convergence.R (the phi_tol checker).
#' @noRd
NULL

#' NA_standards
#'
#' Standards deemed not applicable to PESTO, each with a justification.
#'
#' @srrstatsNA {G1.6} No other R package implements iterative ensemble smoothers, so there is no alternative R implementation to compare against; PESTO is instead compared with the reference external tools PEST 18.25 / pestpp-ies 5.2.16 (see G1.5, G5.4b).
#' @srrstatsNA {G2.4c} PESTO carries no character data requiring conversion via as.character().
#' @srrstatsNA {G2.4d} PESTO accepts no factor inputs, so as.factor() conversion does not apply.
#' @srrstatsNA {G2.4e} PESTO accepts no factor inputs, so conversion from factor does not apply.
#' @srrstatsNA {G2.5} Inputs are numeric matrices / vectors and numeric-column data.frames; factor inputs are not accepted, so ordered-factor handling does not apply.
#' @srrstatsNA {G2.9} PESTO performs no lossy type conversion; the only metadata it adds is default column names when names are absent, which is documented at @param rather than signalled per-call (to avoid noise in ensemble loops).
#' @srrstatsNA {G2.11} PESTO requires numeric parameter columns; a column with a non-standard class carrying a numeric payload (e.g. units) is coerced to its double value during the as.matrix / storage.mode pre-processing.
#' @srrstatsNA {G2.12} List-columns are outside PESTO's numeric-matrix input contract; the numeric coercion rejects them rather than processing them.
#' @srrstatsNA {G2.14c} Imputation is not appropriate for forward-model failures in an ensemble smoother; failed realisations are carried (on_failure = "na") or abort (on_failure = "stop").
#' @srrstatsNA {G5.1} No data/ directory or sysdata.rda; all test data are generated inline from fixed set.seed() calls (fully reproducible from the test sources), so there is no stored data set to export.
#' @srrstatsNA {G5.4c} No canonical published numeric output exists for the iterative ensemble smoother / MDA algorithms to store; correctness is established against an analytic closed-form solution (G5.4) and against fixed-version pestpp-ies 5.2.16 (G5.4b).
#' @srrstatsNA {G5.11} Extended tests use only inline simulated data (cheap, self-contained); no large external data sets or downloads are required.
#' @srrstatsNA {G5.11a} No test performs a download, so there is no download-failure path.
#' @srrstatsNA {BS1.5} A single convergence notion (the phi-reduction rule); there are no multiple convergence checkers to contrast.
#' @srrstatsNA {BS2.9} The IES update holds no internal RNG state and runs no chains -- it is deterministic given the prior ensemble; per-chain seeding does not apply (reproducibility is set by the upstream prior-ensemble draw, tested in test-noise-susceptibility.R).
#' @srrstatsNA {BS2.10} No seed arguments and no computational chains; identical-seed diagnostics do not apply.
#' @srrstatsNA {BS2.14} Warnings flag genuine numerical / specification conditions (e.g. localisation with use_approx = FALSE) and are not silenced by design; suppressWarnings() remains available to the caller.
#' @srrstatsNA {BS3.1} Rank deficiency / collinearity is handled intrinsically by SVD truncation (eigthresh) in the ensemble update; there is no separate collinearity-diagnosis routine.
#' @srrstatsNA {BS3.2} SVD truncation processes (near-)collinear ensembles inherently; no distinct bypass routine is required.
#' @srrstatsNA {BS5.0} The update holds no internal RNG / seed state (deterministic given the prior ensemble); reproducibility is set by the upstream prior draw, so there is no internal seed to return.
#' @srrstatsNA {BS5.2} The assimilation data specification (obs_target, obs_sd, weights) is returned; the prior ensemble and forward model are user-supplied objects retained by the caller and are not copied into the result, avoiding large duplication for expensive simulators.
#' @srrstatsNA {BS5.4} A single convergence checker; there is no checker selection to report.
#' @srrstatsNA {BS6.2} An ensemble smoother has no Markov-chain burn-in; the phi-convergence trace (the default plot() method) is the iteration diagnostic.
#' @srrstatsNA {BS6.4} Optional ("may"); the default print method provides the run summary.
#' @srrstatsNA {BS6.5} Optional ("may").
#' @noRd
NULL
