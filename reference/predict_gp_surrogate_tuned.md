# Predict from an MLE-Tuned GP Surrogate

Companion predictor for
[`train_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/train_gp_surrogate_tuned.md).
It reapplies the per-axis pre-scaling the tuner stored (so an
anisotropic surrogate predicts on the same geometry it was fitted on)
and adds the centred response mean back, returning predictions on the
original response scale.

## Usage

``` r
predict_gp_surrogate_tuned(gp, X_new)
```

## Arguments

- gp:

  A surrogate from
  [`train_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/train_gp_surrogate_tuned.md).

- X_new:

  Numeric matrix of inputs to predict at.

## Value

The list
[`predict_gp_surrogate()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate.md)
returns (`mean`, `variance`, `uncertainty`), with `mean` on the original
response scale.

## See also

[`train_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/train_gp_surrogate_tuned.md).

## Examples

``` r
set.seed(1L)
X <- matrix(runif(40L * 2L), 40L, 2L)
y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2
gp <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))
pred <- predict_gp_surrogate_tuned(gp, X)
```
