# Train a Gaussian Process Surrogate

Trains a GP surrogate model from parameter-observation pairs. Uses
squared exponential (RBF) kernel with automatic relevance determination
(ARD) via median heuristic for length scale.

## Usage

``` r
train_gp_surrogate(
  X_train,
  Y_train,
  length_scale = 0,
  signal_var = 0,
  noise_var = 1e-04
)
```

## Arguments

- X_train:

  Matrix (n x npar). Training parameter sets.

- Y_train:

  Matrix (n x nobs). Corresponding model outputs.

- length_scale:

  Numeric. Kernel length scale. If 0 (default), uses the median
  heuristic (median pairwise distance).

- signal_var:

  Numeric. Signal variance. If 0, uses variance of Y.

- noise_var:

  Numeric. Observation noise variance.

## Value

A list of class `step_gp` containing trained GP components: K_inv
(inverse kernel matrix), alpha (weight vectors), hyperparameters.

## Details

The GP learns the mapping: parameters -\> observations, enabling cheap
prediction of model outputs for new parameter sets.

## Examples

``` r
set.seed(1L)
X_train <- matrix(rnorm(20 * 4), 20, 4)
Y_train <- matrix(rnorm(20 * 6), 20, 6)
gp <- train_gp_surrogate(X_train, Y_train)
pred <- predict_gp_surrogate(gp, X_train)
dim(pred$mean)
#> [1] 20  6
length(pred$uncertainty)
#> [1] 20
```
