# Predict with RFF Sparse GP Surrogate

Predict with RFF Sparse GP Surrogate

## Usage

``` r
predict_rff_surrogate(rff, X_new)
```

## Arguments

- rff:

  A trained RFF model (from `train_rff_surrogate`).

- X_new:

  Matrix (m x npar). New parameter sets.

## Value

A list with mean predictions and approximate uncertainties.

## Examples

``` r
set.seed(1L)
X_train <- matrix(rnorm(30 * 4), 30, 4)
Y_train <- matrix(rnorm(30 * 6), 30, 6)
rff     <- train_rff_surrogate(X_train, Y_train, n_features = 100L)
X_new   <- matrix(rnorm(5 * 4), 5, 4)
pred    <- predict_rff_surrogate(rff, X_new)
dim(pred$mean)
#> [1] 5 6
length(pred$uncertainty)
#> [1] 5
```
