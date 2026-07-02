# Ensemble Solution with Adaptive SVD Backend

A variant of
[`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
that selects the SVD backend automatically – a randomised SVD for
low-rank problems, otherwise a dense LAPACK / Accelerate decomposition –
and returns the upgrade together with timing diagnostics. It is the
convenient entry point when the best backend for a given problem size is
not known in advance. All computation is on the CPU.

## Usage

``` r
ensemble_solution_adaptive(
  par_diff,
  obs_diff,
  obs_resid,
  par_resid,
  weights,
  parcov_inv,
  Am,
  cur_lam,
  eigthresh = 1e-06,
  use_approx = TRUE,
  use_prior_scaling = FALSE,
  iter = 1L,
  reg_factor = -1,
  svd_method = "auto",
  target_rank = 0L
)
```

## Arguments

- par_diff:

  Matrix (npar x nreal). Parameter anomalies.

- obs_diff:

  Matrix (nobs x nreal). Observation anomalies.

- obs_resid:

  Matrix (nobs x nreal). Observation residuals.

- par_resid:

  Matrix (npar x nreal). Parameter residuals.

- weights:

  Numeric vector (nobs). Observation weights.

- parcov_inv:

  Numeric vector (npar). Inverse parameter covariance diagonal.

- Am:

  Matrix. Random Am matrix for upgrade_2.

- cur_lam:

  Numeric. Marquardt lambda.

- eigthresh:

  Numeric. Eigenvalue truncation threshold.

- use_approx:

  Logical. Skip upgrade_2.

- use_prior_scaling:

  Logical. Scale by prior covariance.

- iter:

  Integer. Current iteration.

- reg_factor:

  Numeric. Regularisation factor.

- svd_method:

  Character. SVD method: "auto", "rsvd", "accelerate", "eigen".

- target_rank:

  Integer. Target rank for randomised SVD (0 = auto).

## Value

A list with upgrade matrix and performance diagnostics.

## Examples

``` r
set.seed(1L)
npar  <- 4L
nreal <- 20L
nobs  <- 30L
par_diff  <- matrix(rnorm(npar * nreal), npar, nreal)
obs_diff  <- matrix(rnorm(nobs * nreal), nobs, nreal)
obs_resid <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
par_resid <- matrix(rnorm(npar * nreal, sd = 0.1), npar, nreal)
weights    <- rep(1, nobs)
parcov_inv <- rep(1, npar)
Am         <- matrix(0, 0, 0)
res <- ensemble_solution_adaptive(
  par_diff   = par_diff,
  obs_diff   = obs_diff,
  obs_resid  = obs_resid,
  par_resid  = par_resid,
  weights    = weights,
  parcov_inv = parcov_inv,
  Am         = Am,
  cur_lam    = 1.0,
  svd_method = "auto"
)
dim(res$upgrade)
#> [1] 20  4
```
