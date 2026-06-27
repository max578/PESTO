#' srr_stats
#'
#' Compliance documentation for the rOpenSci statistical-software standards.
#' `@srrstats` tags lead with how PESTO complies; `@srrstatsTODO` marks standards
#' not yet addressed; non-applicable standards are in the NA_standards block.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} Primary references are cited in roxygen (Chen & Oliver 2013; Evensen 2018; Goudriaan & Monteith 1990; Anderson & May 1991).
#' @srrstatsTODO {G1.1} *Statistical Software should document whether the algorithm(s) it implements are:* - *The first implementation of a novel algorithm*; or - *The first implementation within **R** of an algorithm which has previously been implemented in other languages or contexts*; or - *An improvement on other implementations of similar algorithms in **R***.
#' @srrstats {G1.2} A Life Cycle Statement is provided in API_STABILITY.md (additive-then-lifecycle policy).
#' @srrstatsTODO {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G1.4} All exported functions are documented with roxygen2 (see man/).
#' @srrstats {G1.4a} Internal (non-exported) helpers carry @noRd and @keywords internal (e.g. R/internal_validation.R).
#' @srrstatsTODO {G1.5} *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.*
#' @srrstatsTODO {G1.6} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.*
#' @srrstats {G2.0} Input lengths are asserted via the .assert_*_scalar helpers (length(x) != 1L) in R/internal_validation.R.
#' @srrstatsTODO {G2.0a} Provide explicit secondary documentation of any expectations on lengths of inputs
#' @srrstats {G2.1} Input types are asserted (numeric / character / logical / function / matrix) in R/internal_validation.R.
#' @srrstatsTODO {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#' @srrstats {G2.2} Univariate parameters are guarded by scalar length assertions (R/internal_validation.R).
#' @srrstats {G2.3} Univariate character inputs are restricted via match.arg() / .assert_choice().
#' @srrstats {G2.3a} match.arg() is used for character-parameter dispatch (method, format, on_failure, taper).
#' @srrstats {G2.3b} Unmatched character values raise an explicit error via match.arg().
#' @srrstatsTODO {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstatsTODO {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstatsTODO {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstatsTODO {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstatsTODO {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstatsTODO {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstatsTODO {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.*
#' @srrstatsTODO {G2.7} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.*
#' @srrstatsTODO {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#' @srrstatsTODO {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost or added.*
#' @srrstatsTODO {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs behaves consistently regardless of the class of tabular data used as input.*
#' @srrstatsTODO {G2.11} *Software should ensure that data.frame-like tabular objects with columns lacking standard class attributes are appropriately processed, and this behaviour should be tested.*
#' @srrstatsTODO {G2.12} *Software should ensure that data.frame-like tabular objects with list columns are appropriately pre-processed, and this behaviour should be tested.*
#' @srrstatsTODO {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstatsTODO {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#' @srrstatsTODO {G2.14a} *error on missing data*
#' @srrstatsTODO {G2.14b} *ignore missing data with default warnings or messages issued*
#' @srrstatsTODO {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstatsTODO {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default na.rm = FALSE-type parameters.*
#' @srrstatsTODO {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
#' @srrstats {G3.0} No naive floating-point equality is used (verified by source scan); tolerances / inequalities are used throughout.
#' @srrstatsTODO {G3.1} *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
#' @srrstatsTODO {G3.1a} *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).*
#' @srrstatsTODO {G4.0} *Statistical Software which enables outputs to be written to local files should parse parameters specifying file names to ensure appropriate file suffixes are automatically generated where not provided.*
#'
#' Testing standards (the G5 series) are documented next to the tests that
#' satisfy them, under tests/testthat/. Non-applicable members (G5.1, G5.4c,
#' G5.11, G5.11a) are in the NA_standards block below.
#'
#' Bayesian and Monte Carlo (BS) standards -- compliance leads each entry.
#' @srrstats {BS1.0} "hyperparameter" is used only for the GP-surrogate kernel parameters (length-scale, input-scale), defined in train_gp_surrogate_tuned() documentation.
#' @srrstats {BS1.1} Data entry (a prior parameter ensemble, a named observation vector, obs_sd) is shown in the pesto_ies_callback() examples and the getting-started vignette.
#' @srrstats {BS1.2} The prior is supplied as a prior parameter ensemble; how to specify it is documented at the prior_ensemble parameter with worked examples.
#' @srrstatsTODO {BS1.2a} Prior specification in the main README (pending the README.Rmd conversion).
#' @srrstats {BS1.2b} The getting-started vignette shows prior vs posterior ensembles via plot_ensemble().
#' @srrstats {BS1.2c} pesto_ies_callback() examples construct a prior ensemble and run the smoother.
#' @srrstats {BS1.3} The computational-process parameters (noptmax, lambda, parcov, eigthresh, use_approx) are each documented in roxygen.
#' @srrstatsTODO {BS1.3a} Document (text + example) using a previous run's posterior as the starting ensemble of the next run.
#' @srrstats {BS1.3b} The alternative drivers (in-process callback, sequential filter, MDA, .pst path) are documented and cross-linked via @seealso.
#' @srrstats {BS2.1} Dimensional commensurability of inputs (nobs / npar / nreal) is enforced in the driver pre-processing.
#' @srrstats {BS2.1a} The dimension / length validation is tested (test-ies-callback.R, test-edge-conditions.R).
#' @srrstats {BS2.2} obs_sd (the observation-error specification) is validated as a distinct pre-processing step before the iteration loop.
#' @srrstats {BS2.3} obs_sd must be scalar or length-nobs; a non-conforming length raises an error (no excess values silently discarded).
#' @srrstats {BS2.4} obs_sd length is checked for commensurability with the observations (scalar recycled to nobs).
#' @srrstats {BS2.5} Variance-like parameters are constrained positive: obs_sd > 0 and parcov > 0 are enforced.
#' @srrstats {BS2.6} Computational parameters are range-checked (noptmax >= 1, positive eigthresh, in-range fidelity levels).
#' @srrstats {BS2.11} Starting points are supplied as the prior_ensemble matrix -- inherently plural, one row per realisation.
#' @srrstats {BS2.12} The verbose parameter controls progress output and defaults to TRUE.
#' @srrstats {BS2.15} on_failure = "na" catches forward-model failures and captures them in the return value (failure_rate + NA realisations) rather than aborting; tested in test-ies-callback.R.
#' @srrstats {BS3.0} Missing-value handling is documented via the on_failure parameter and the forward-model contract (failed realisations become NA rows, tolerated or aborted).
#' @srrstats {BS4.0} The assimilation algorithms are documented via citation (Chen & Oliver 2013; Evensen 2018) in roxygen references.
#' @srrstatsTODO {BS4.3} Offer a convergence checker (e.g. a phi-reduction stopping rule) with a documented reference -- design decision pending (fixed-iteration vs adaptive stopping).
#' @srrstatsTODO {BS4.4} Enable stopping on convergence -- coupled to the BS4.3 design decision.
#' @srrstats {BS4.5} Non-convergence is surfaced via the returned phi trace and the per-iteration spread-ESS ratio; failed forward evaluations are handled via on_failure; Marquardt lambda damping mitigates divergence.
#'
#' Bayesian standards documented next to their code/tests: BS2.7, BS2.8, BS2.13,
#' BS5.1, BS5.3, BS5.5 (test-bayesian-interface.R); BS4.1, BS7.2
#' (test-correctness-analytic.R); BS4.2, BS7.0, BS7.1, BS7.4, BS7.4a
#' (test-bayesian-recovery.R); BS7.3 (test-recovery-scaling.R); BS6.0, BS6.1,
#' BS6.3 (R/ies_result_methods.R).
#' @noRd
NULL

#' NA_standards
#'
#' Standards deemed not applicable to PESTO, each with a justification.
#'
#' @srrstatsNA {G2.5} PESTO inputs are numeric matrices / vectors and numeric-column data.frames; factor inputs are not accepted, so ordered-factor handling does not apply.
#' @srrstatsNA {G5.1} No data/ directory or sysdata.rda; all test data are generated inline from fixed set.seed() calls (fully reproducible from the test sources), so there is no stored data set to export.
#' @srrstatsNA {G5.4c} No canonical published numeric output exists for the iterative ensemble smoother / MDA algorithms to store; correctness is established against an analytic closed-form solution (G5.4) and against fixed-version pestpp-ies 5.2.16 (G5.4b).
#' @srrstatsNA {G5.11} Extended tests use only inline simulated data (cheap, self-contained); no large external data sets or downloads are required.
#' @srrstatsNA {G5.11a} No test performs a download, so there is no download-failure path.
#' @srrstatsNA {BS1.4} The iterative ensemble smoother runs a fixed number of iterations with no in-run convergence checker; documenting use with/without a checker does not apply.
#' @srrstatsNA {BS1.5} A single convergence notion (the phi misfit); there are no multiple convergence checkers to contrast.
#' @srrstatsNA {BS2.9} The IES update holds no internal RNG state and runs no chains -- it is deterministic given the prior ensemble; per-chain seeding does not apply (reproducibility is set by the upstream prior-ensemble draw, tested in test-noise-susceptibility.R).
#' @srrstatsNA {BS2.10} No seed arguments and no computational chains; identical-seed diagnostics do not apply.
#' @srrstatsNA {BS2.14} Warnings flag genuine numerical / specification conditions (e.g. localisation with use_approx = FALSE) and are not silenced by design; suppressWarnings() remains available to the caller.
#' @srrstatsNA {BS3.1} Rank deficiency / collinearity is handled intrinsically by SVD truncation (eigthresh) in the ensemble update; there is no separate collinearity-diagnosis routine.
#' @srrstatsNA {BS3.2} SVD truncation processes (near-)collinear ensembles inherently; no distinct bypass routine is required.
#' @srrstatsNA {BS4.6} No in-run convergence checker, so checker-vs-fixed equivalence cannot be (and need not be) tested.
#' @srrstatsNA {BS4.7} No parametrised convergence checker whose threshold effects could be tested.
#' @srrstatsNA {BS5.0} The update holds no internal RNG / seed state (deterministic given the prior ensemble); reproducibility is set by the upstream prior draw, so there is no internal seed to return.
#' @srrstatsNA {BS5.2} The assimilation data specification (obs_target, obs_sd, weights) is returned; the prior ensemble and forward model are user-supplied objects retained by the caller and are not copied into the result, avoiding large duplication for expensive simulators.
#' @srrstatsNA {BS5.4} A single convergence notion; there is no checker selection to report.
#' @srrstatsNA {BS6.2} An ensemble smoother has no Markov-chain burn-in; the phi-convergence trace (the default plot() method) is the iteration diagnostic.
#' @srrstatsNA {BS6.4} Optional ("may"); the default print method provides the run summary.
#' @srrstatsNA {BS6.5} Optional ("may").
#' @noRd
NULL
