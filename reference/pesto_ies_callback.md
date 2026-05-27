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
  parcov = NULL,
  eigthresh = 1e-06,
  use_approx = TRUE,
  on_failure = c("na", "stop"),
  verbose = TRUE
)
```

## Arguments

- forward_model:

  A function with signature `function(theta) -> obs`, where `theta` is
  an `nreal x npar` numeric matrix and `obs` is an `nreal x nobs`
  numeric matrix. Failed realisations may return rows of `NA`; the
  driver tolerates them (see `on_failure`).

- prior_ensemble:

  Matrix or data.table, `nreal x npar`. Columns are parameters; an
  optional `real_name` column is preserved if present. Column names
  supply parameter names.

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

- parcov:

  Numeric vector of length `npar`, the diagonal of the prior parameter
  covariance. Defaults to the column-wise variance of `prior_ensemble`;
  zero or negative entries are replaced with 1.0.

- eigthresh:

  Numeric. SVD eigenvalue truncation threshold (default 1e-6).

- use_approx:

  Logical. If TRUE (default), skip the prior-scaling correction
  (upgrade_2); matches the typical `pestpp-ies` default.

- on_failure:

  Character. `"na"` (default) carries failed realisations forward
  unchanged and proceeds; `"stop"` aborts on any failure.

- verbose:

  Logical. Print per-iteration phi summaries.

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

  List of per-iteration metadata (lambda, mean phi, failure count).

- runtime_seconds:

  Total wall-clock runtime.

- n_forward_evals:

  Total number of realisation-level forward evaluations across all
  iterations (including the final refresh).

- failure_rate:

  Fraction of forward evaluations that returned NA.

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

## References

Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
iterative ensemble smoother for efficient history matching and
uncertainty quantification. *Computational Geosciences*, 17(4), 689–703.

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
