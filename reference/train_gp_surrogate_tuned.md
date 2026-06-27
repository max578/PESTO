# Train a GP Surrogate with Maximum-Likelihood (Anisotropic) Length Scales

[`train_gp_surrogate()`](https://max578.github.io/PESTO/reference/train_gp_surrogate.md)
defaults to a single median-heuristic length scale: fast and robust, but
on a strongly anisotropic or multi-scale response it can be several
times less accurate than length scales tuned to the data. This helper
keeps the same fast C++ GP but selects the length scale(s) by
**maximising the GP's own log marginal likelihood** – the criterion a
maximum-likelihood Gaussian-process library (for example `DiceKriging`)
optimises.

## Usage

``` r
train_gp_surrogate_tuned(
  X_train,
  Y_train,
  anisotropic = TRUE,
  signal_var = NULL,
  noise_var = 1e-04,
  n_restarts = 5L,
  n_grid = 40L,
  length_scale_bounds = NULL
)
```

## Arguments

- X_train:

  Numeric matrix of training inputs (rows = points, columns =
  parameters).

- Y_train:

  Numeric matrix (or vector) of training outputs.

- anisotropic:

  Logical. If `TRUE` (default) and there are at least two input
  dimensions, estimate one length scale per dimension; otherwise a
  single length scale.

- signal_var:

  Numeric or `NULL` (default, the GP's automatic mean per-output
  response variance).

- noise_var:

  Numeric observation-noise variance / nugget. Default `1e-4`.

- n_restarts:

  Integer. Random restarts for the anisotropic optimisation; the best
  marginal likelihood is kept. Default `5`.

- n_grid:

  Integer. Length scales in the isotropic grid search. Default `40`.

- length_scale_bounds:

  Numeric `c(lower, upper)` for the isotropic search, or `NULL`
  (default) to derive it from the pairwise distances of `X_train`.

## Value

The list from
[`train_gp_surrogate()`](https://max578.github.io/PESTO/reference/train_gp_surrogate.md)
at the maximum-likelihood length scale(s), trained on centred (and, when
anisotropic, pre-scaled) data, with an added `tuning` element:
`anisotropic`, `length_scale` (per dimension when anisotropic, scalar
otherwise), `input_scale` (the per-axis divisor applied before training;
all ones when isotropic), `y_mean` (the per-output centring offset),
`length_scale_median`, `log_marginal_likelihood` (at the optimum) and
`log_marginal_likelihood_median` (at the single median heuristic, the
baseline this improves on).

## Details

By default the fit is **anisotropic**: one length scale is estimated per
input dimension by pre-scaling each coordinate (the isotropic C++ kernel
applied to coordinates divided by per-axis length scales is an
anisotropic kernel on the originals), optimised with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) from several
starts. On a strongly anisotropic response this recovers most of the
accuracy a single length scale leaves on the table: on the Branin
function it cuts held-out error roughly three-fold and brings the
surrogate to within a small factor of an anisotropic MLE oracle, where a
single length scale sits about seven-fold worse. With
`anisotropic = FALSE`, or a single input dimension, one length scale is
tuned by a log-spaced grid plus a
[`stats::optimize()`](https://rdrr.io/r/stats/optimize.html) refinement.

The response is centred before fitting (the C++ GP is zero-mean), and
the marginal variance is left at the GP's automatic value unless
`signal_var` is supplied.

Because an anisotropic fit stores the GP on pre-scaled coordinates,
**predict with
[`predict_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate_tuned.md)**,
which reapplies the pre-scaling and adds the centred mean back. Calling
[`predict_gp_surrogate()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate.md)
directly on a tuned surrogate is correct only in the isotropic case.

## See also

[`predict_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate_tuned.md)
to predict from the result;
[`train_gp_surrogate()`](https://max578.github.io/PESTO/reference/train_gp_surrogate.md)
for the default single-heuristic fit.

## Examples

``` r
set.seed(1L)
X <- matrix(runif(40L * 2L), 40L, 2L)
y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2
gp <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))
gp$tuning$length_scale
#> [1] 0.3634694 0.9101072
pred <- predict_gp_surrogate_tuned(gp, X)
```
