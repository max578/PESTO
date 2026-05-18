# Ensemble Solution Kernel (GLM form)

Implements the core IES ensemble update equation using the
Gauss-Levenberg-Marquardt (GLM) formulation from Chen & Oliver (2013).
This is the computational hotspot of the iterative ensemble smoother.

## Usage

``` r
ensemble_solution(
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
  reg_factor = -1
)
```

## Arguments

- par_diff:

  Matrix (npar x nreal). Parameter anomalies (deviations from mean).

- obs_diff:

  Matrix (nobs x nreal). Observation anomalies.

- obs_resid:

  Matrix (nobs x nreal). Observation residuals, simulated minus observed
  (sim - obs). This sign convention is required so that the leading
  negative in the GLM update (\\\Delta\theta = -\Delta\theta' V s (s^2 +
  \lambda I)^{-1} U^T r\\) yields a descent step on \\\Phi =
  \\W(g(\theta) - d^{obs})\\^2\\. Passing (obs - sim) inverts the
  gradient and causes phi to diverge.

- par_resid:

  Matrix (npar x nreal). Parameter residuals (par - prior mean).

- weights:

  Numeric vector (nobs). Observation weights (1/sqrt(variance)).

- parcov_inv:

  Numeric vector (npar). Diagonal of inverse parameter covariance.

- Am:

  Matrix (npar x nreal-1). Random Am matrix for upgrade_2 (optional).

- cur_lam:

  Numeric. Current Marquardt lambda.

- eigthresh:

  Numeric. Eigenvalue truncation threshold (0-1).

- use_approx:

  Logical. If TRUE, skip upgrade_2 (prior-scaling correction).

- use_prior_scaling:

  Logical. Scale by prior covariance.

- iter:

  Integer. Current iteration number.

- reg_factor:

  Numeric. Regularisation factor for upgrade_2 blending.

## Value

Matrix (nreal x npar). Parameter upgrade vectors (one row per
realisation). The returned matrix is the *negative-direction* step from
the Chen-Oliver 2013 GLM update formula. To advance the ensemble, apply
by subtraction: `par_new = par_old - upgrade`. The R-side driver
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
handles this convention internally.

## Details

The update equation solves: \$\$\Delta\theta = -\Delta\theta' V s (s^2 +
\lambda I)^{-1} U^T r\$\$ where the SVD is performed on the scaled
observation difference matrix.

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
upgrade <- ensemble_solution(
  par_diff   = par_diff,
  obs_diff   = obs_diff,
  obs_resid  = obs_resid,
  par_resid  = par_resid,
  weights    = weights,
  parcov_inv = parcov_inv,
  Am         = Am,
  cur_lam    = 1.0
)
dim(upgrade)
#> [1] 20  4
```
