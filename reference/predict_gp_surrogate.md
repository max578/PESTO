# Predict with GP Surrogate (with Uncertainty)

Generates predictions and prediction uncertainties for new parameter
sets using a trained GP surrogate. The uncertainty estimates are crucial
for the adaptive switching criterion.

## Usage

``` r
predict_gp_surrogate(gp, X_new)
```

## Arguments

- gp:

  A trained GP model (from `train_gp_surrogate`).

- X_new:

  Matrix (m x npar). New parameter sets to predict.

## Value

A list with:

- mean:

  Matrix (m x nobs). Predicted observations.

- variance:

  Matrix (m x nobs). Prediction variance per output.

- uncertainty:

  Numeric vector (m). Mean prediction uncertainty per realisation.

## Examples

``` r
set.seed(1L)
X_train <- matrix(rnorm(20 * 4), 20, 4)
Y_train <- matrix(rnorm(20 * 6), 20, 6)
gp      <- train_gp_surrogate(X_train, Y_train)
X_new   <- matrix(rnorm(5 * 4), 5, 4)
pred    <- predict_gp_surrogate(gp, X_new)
dim(pred$mean)
#> [1] 5 6
length(pred$uncertainty)
#> [1] 5
```
