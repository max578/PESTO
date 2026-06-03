# Run a Sequential (Filter-Mode) Iterative Ensemble Smoother

Assimilates observations in time-ordered **windows** against a static
parameter ensemble. Each window runs an IES Gauss-Levenberg-Marquardt
update using only that window's observation block; the updated ensemble
is carried forward as the prior for the next window. This is the
filtering analogue of
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md):
instead of one batch assimilation of all observations, the posterior is
refined window by window and is available after each.

## Usage

``` r
pesto_ies_filter(
  forward_model,
  prior_ensemble,
  obs,
  obs_sd,
  windows,
  window_noptmax = 1L,
  lambda = 1,
  fidelity_schedule = NULL,
  parcov = NULL,
  eigthresh = 1e-06,
  use_approx = TRUE,
  inflation = NULL,
  localisation = NULL,
  on_failure = c("na", "stop"),
  verbose = TRUE
)
```

## Arguments

- forward_model:

  A function `function(theta) -> obs`, a
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md),
  or a
  [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  (then see `fidelity_schedule`). The model maps the full `nreal x npar`
  parameter matrix to the full `nreal x nobs` observation matrix;
  windows select columns of that output. A bare function is auto-wrapped
  via
  [`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md).
  See
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  for the contract.

- prior_ensemble:

  Matrix or data.table, `nreal x npar`.

- obs:

  Named numeric vector of length `nobs`. The full set of target
  observations; `windows` indexes into it.

- obs_sd:

  Numeric scalar or length-`nobs` vector. Observation standard
  deviation(s); weights are `1 / obs_sd`.

- windows:

  A list of integer vectors. Each element gives the observation indices
  (into `obs`, `1`-based) assimilated at that window, in assimilation
  order. Windows must be **disjoint** (no observation assimilated twice)
  and need not cover all of `obs`.

- window_noptmax:

  Integer scalar or vector. IES iterations *within* each window (default
  `1L`, the pure filter; `> 1` gives an iterated filter per window). A
  scalar is recycled; a short vector is right-padded with its last
  value.

- lambda:

  Numeric scalar or per-window vector. Marquardt lambda (default `1.0`;
  recycled / right-padded across windows).

- fidelity_schedule:

  Integer vector or `NULL`. Per-window fidelity level for a
  [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  (see *Multi-fidelity*).

- parcov:

  Numeric vector of length `npar`, the diagonal of the prior parameter
  covariance. Defaults to the prior ensemble's column-wise variance
  (non-positive entries replaced with `1`).

- eigthresh:

  Numeric. SVD eigenvalue truncation (default `1e-6`).

- use_approx:

  Logical. If `TRUE` (default) skip the prior-scaling correction,
  matching the `pestpp-ies` default.

- inflation:

  A
  [`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)
  specification, or `NULL` (default). Applied within each window's inner
  update to counteract ensemble under-dispersion. See
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).

- localisation:

  A
  [`pesto_localisation()`](https://max578.github.io/PESTO/reference/pesto_localisation.md)
  specification, or `NULL` (default). Tapers the per-window Kalman gain
  to suppress spurious finite-sample correlations. See
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).

- on_failure:

  Character. `"na"` (default) tolerates failed realisations; `"stop"`
  aborts on any failure.

- verbose:

  Logical. Print a per-window phi summary.

## Value

A list of class `c("pesto_ies_filter_result", "pesto_ies_result")` with
components:

- phi:

  data.table of per-realisation phi by window (on that window's
  observation block).

- par_ensemble:

  Final parameter ensemble (data.table).

- obs_ensemble:

  Final simulated-observation ensemble, full `nobs` columns
  (data.table).

- windows:

  List of per-window metadata: assimilated indices, lambda, mean phi,
  per-parameter ensemble mean and standard deviation (the sd trace shows
  the posterior tightening), and failure count.

- runtime_seconds, n_forward_evals, failure_rate:

  Run totals.

- fidelity:

  Multi-fidelity provenance (or `NULL`), as in
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).

## Filter vs smoother

[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
forms parameter residuals against the *original* prior mean and
assimilates every observation together. Here each window forms residuals
against the *current* (carried-forward) ensemble mean and assimilates
only its own block, so information accrues sequentially. With the
default `use_approx = TRUE` the carried-forward ensemble itself encodes
the background covariance, so the sequential behaviour comes from
propagating the ensemble, not from an explicit prior term.

## Multi-fidelity

When `forward_model` is a
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md),
`fidelity_schedule` selects the fidelity level evaluated at each
*window* (recycled / padded to the number of windows; default: highest
fidelity throughout). The final ensemble refresh always uses the highest
fidelity.

## References

Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
iterative ensemble smoother for efficient history matching and
uncertainty quantification. *Computational Geosciences*, 17(4), 689–703.

## See also

[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
for the batch smoother;
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
for fidelity stacks;
[`as_manifest()`](https://max578.github.io/PESTO/reference/as_manifest.md)
to wrap the result in the ensemble-manifest contract.

## Examples

``` r
# Linear-Gaussian recovery, assimilated in three observation windows.
set.seed(1)
npar <- 3; nobs <- 9; nreal <- 80
G <- matrix(rnorm(nobs * npar), nobs, npar)
theta_true <- c(1.0, -0.5, 2.0)
y <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = 0.05)
f <- function(theta) theta %*% t(G)
prior <- matrix(rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, paste0("p", 1:npar)))
fit <- pesto_ies_filter(
  forward_model  = f, prior_ensemble = prior,
  obs = setNames(y, paste0("o", 1:nobs)), obs_sd = 0.05,
  windows = list(1:3, 4:6, 7:9), verbose = FALSE
)
# Posterior sd should shrink window over window:
vapply(fit$windows, function(w) mean(w$par_sd), numeric(1))
#> [1] 0.052603202 0.003631259 0.003615228
```
