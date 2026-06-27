# Run IES with an In-Process R Callback Forward Model

Drives an Iterative Ensemble Smoother entirely in R, using a
user-supplied forward model callable instead of the PEST++ `.pst`-file
invocation cycle. Each iteration:

## Usage

``` r
pesto_ies_callback(
  forward_model,
  prior_ensemble,
  obs,
  obs_sd,
  noptmax = 4L,
  lambda = 1,
  fidelity_schedule = NULL,
  parcov = NULL,
  eigthresh = 1e-06,
  use_approx = TRUE,
  inflation = NULL,
  localisation = NULL,
  on_failure = c("na", "stop"),
  verbose = TRUE,
  phi_tol = NULL
)
```

## Arguments

- forward_model:

  One of: a function with signature `function(theta) -> obs` (where
  `theta` is an `nreal x npar` numeric matrix and `obs` an
  `nreal x nobs` numeric matrix); a
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  (the typed contract object, e.g. one carrying a
  `parallel = "multicore"` evaluation strategy); or a
  [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  (then see `fidelity_schedule`). A bare function is auto-wrapped via
  [`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md)
  with this call's `on_failure`. Failed realisations may return rows of
  `NA`; the driver tolerates them (see `on_failure`). To run
  realisations in parallel, pass a `pesto_forward_model` built with
  `parallel = "multicore"` rather than a bare function.

- prior_ensemble:

  Matrix or data.table, `nreal x npar`. Columns are parameters; an
  optional `real_name` column is preserved if present. Column names
  supply parameter names. To **warm-start**, pass the parameter columns
  of a previous run's `par_ensemble` as the prior of the next: a
  posterior ensemble is itself a valid prior ensemble.

- obs:

  Named numeric vector. Target observations.

- obs_sd:

  Numeric scalar or vector of length `nobs`. Observation standard
  deviation(s); the IES weights are `1/obs_sd`.

- noptmax:

  Integer. Number of IES iterations (default 4).

- lambda:

  Numeric scalar or vector. Marquardt lambda per iteration. A scalar is
  recycled; a vector shorter than `noptmax` is right-padded with its
  last value (default 1.0).

- fidelity_schedule:

  Integer vector or `NULL`. Only consulted when `forward_model` is a
  [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md):
  the fidelity level to evaluate at each iteration (recycled /
  right-padded to `noptmax`). `NULL` (default) uses the highest fidelity
  every iteration. Ignored otherwise.

- parcov:

  Numeric vector of length `npar`, the diagonal of the prior parameter
  covariance. Defaults to the column-wise variance of `prior_ensemble`;
  zero or negative entries are replaced with 1.0.

- eigthresh:

  Numeric. SVD eigenvalue truncation threshold (default 1e-6).

- use_approx:

  Logical. If TRUE (default), skip the prior-scaling correction
  (upgrade_2); matches the typical `pestpp-ies` default.

- inflation:

  A
  [`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)
  specification, or `NULL` (default) for no covariance inflation.
  Inflation counteracts ensemble under-dispersion – the progressive
  collapse of posterior spread a finite-ensemble smoother suffers.
  `NULL` leaves the update identical to the un-inflated smoother.

- localisation:

  A
  [`pesto_localisation()`](https://max578.github.io/PESTO/reference/pesto_localisation.md)
  specification, or `NULL` (default) for no localisation. Localisation
  tapers the Kalman gain to suppress spurious finite-sample
  parameter-observation correlations; the active path uses
  [`ensemble_solution_localised()`](https://max578.github.io/PESTO/reference/ensemble_solution_localised.md)
  (approximate / upgrade_1 form), so a non-`NULL` `localisation` with
  `use_approx = FALSE` warns and drops the null-space correction.

- on_failure:

  Character. `"na"` (default) carries failed realisations forward
  unchanged and proceeds; `"stop"` aborts on any failure.

- verbose:

  Logical. Print per-iteration phi summaries.

- phi_tol:

  Numeric scalar or `NULL`. Optional convergence tolerance: when
  non-`NULL`, iteration stops early once the relative reduction in the
  mean objective function (phi) between successive iterations falls
  below `phi_tol` – the phi-reduction stopping rule of White (2018).
  `NULL` (default) runs the full `noptmax` iterations, leaving the
  update byte-identical to the unchecked smoother. Use a smaller
  `phi_tol` to demand more iterations, a larger one to stop sooner.

## Value

A list of class `c("pesto_ies_callback_result", "pesto_ies_result")`
with components:

- phi:

  data.table of per-realisation phi by iteration.

- par_ensemble:

  Final parameter ensemble (data.table).

- obs_ensemble:

  Final simulated-observation ensemble (data.table).

- iterations:

  List of per-iteration metadata: `lambda`, `mean_phi`, `n_failures`,
  and the dispersion diagnostics `spread_ess` / `spread_ess_ratio`
  ([`ensemble_spread_ess()`](https://max578.github.io/PESTO/reference/ensemble_spread_ess.md)),
  plus `inflation_method` / `inflation_factor` / `retention` and
  `localisation` / `loc_threshold` / `loc_frac_active` when those
  countermeasures are active.

- runtime_seconds:

  Total wall-clock runtime.

- n_forward_evals:

  Total number of realisation-level forward evaluations across all
  iterations (including the final refresh).

- failure_rate:

  Fraction of forward evaluations that returned NA.

- converged:

  Logical: `TRUE` if the run stopped early on the `phi_tol` convergence
  criterion; `FALSE` when `phi_tol` is `NULL` or the full `noptmax` was
  reached.

- n_iterations:

  Number of IES iterations actually run (fewer than `noptmax` if the
  convergence checker stopped the run early).

- fidelity:

  For a
  [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  run, a provenance list
  `list(type, schedule, final_level, n_levels, costs)` recording the
  realised per-iteration fidelity schedule; `NULL` for a single-fidelity
  run. Consumed by
  [`as_manifest()`](https://max578.github.io/PESTO/reference/as_manifest.md)
  to populate the manifest contract.

## Details

1.  Evaluates `forward_model(par_ensemble)` to obtain simulated
    observations.

2.  Forms parameter / observation anomalies and residuals.

3.  Calls the C++ kernel
    [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
    (GLM form, Chen & Oliver 2013) to compute an `nreal x npar` upgrade.

4.  Adds the upgrade to the current ensemble.

The classic `.pst`-file path remains available via
[`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md)
for full PEST++ compatibility. Use this callback driver when the forward
model is itself an R function (e.g. an `apsimx` wrapper from
[`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md),
a Python bridge, or a synthetic test problem) and the per-realisation
file-I/O overhead of the `.pst` path is the bottleneck.

Phase-1 behaviour: single lambda per iteration (or a user-supplied
schedule). A line-search over `ies_lambda_mults` matching `pestpp-ies`
is a planned Phase-2 enhancement; for the common case of a well-behaved
forward model with `lambda = 1`, the GLM update reduces phi reliably
(see vignette `apsim-callback`).

## Multi-fidelity

When `forward_model` is a
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md),
`fidelity_schedule` selects which fidelity level each iteration
evaluates – a recycled / padded integer vector, one entry per iteration
(default: the highest fidelity every iteration, i.e. exactly the
single-fidelity behaviour). This supports fidelity ramping (cheap early
iterations, expensive late ones); the final ensemble refresh always uses
the highest fidelity so the returned posterior is at full resolution.
The control-variate combiner
[`mf_control_variate()`](https://max578.github.io/PESTO/reference/mf_control_variate.md)
is the plug-in point for bias-corrected surrogate cascades.

## References

Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
iterative ensemble smoother for efficient history matching and
uncertainty quantification. *Computational Geosciences*, 17(4), 689–703.

White, J.T. (2018). A model-independent iterative ensemble smoother for
efficient history-matching and uncertainty quantification in very high
dimensions. *Environmental Modelling & Software*, 109, 191–201.

## See also

[`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md)
for the `.pst`-file path;
[`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
for the apsimx adapter.

## Examples

``` r
# Linear-Gaussian recovery toy
set.seed(1)
npar <- 3; nobs <- 6; nreal <- 80
G <- matrix(rnorm(nobs * npar), nobs, npar)
theta_true <- c(1.0, -0.5, 2.0)
y <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = 0.05)
f <- function(theta) theta %*% t(G)
prior <- matrix(rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, paste0("p", 1:npar)))
fit <- pesto_ies_callback(
  forward_model = f, prior_ensemble = prior,
  obs = setNames(y, paste0("o", 1:nobs)), obs_sd = 0.05,
  noptmax = 5, verbose = FALSE
)
colMeans(as.matrix(fit$par_ensemble[, -1]))  # should approach theta_true
#>         p1         p2         p3 
#>  1.0136905 -0.4936753  1.9847267 
```
