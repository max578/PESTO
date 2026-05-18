# Train a Sparse GP Surrogate via Random Fourier Features

Approximates the RBF kernel GP using random Fourier features, reducing
training cost from O(n^3) to O(n \* D^2) where D is the number of random
features (typically 100-500). This enables GP surrogates for ensembles
of 1,000+ realisations.

## Usage

``` r
train_rff_surrogate(
  X_train,
  Y_train,
  n_features = 200L,
  length_scale = 0,
  noise_var = 1e-04
)
```

## Arguments

- X_train:

  Matrix (n x npar). Training parameter sets.

- Y_train:

  Matrix (n x nobs). Corresponding model outputs.

- n_features:

  Integer. Number of random Fourier features (default 200).

- length_scale:

  Numeric. Kernel length scale (0 = median heuristic).

- noise_var:

  Numeric. Observation noise variance.

## Value

A list containing the trained RFF model.

## Details

The RBF kernel k(x,x') = sigma^2 exp(-\|\|x-x'\|\|^2 / 2l^2) is
approximated by k(x,x') ~ z(x)^T z(x') where z(x) = sqrt(2/D) \*
cos(W*x + b), where W are random frequencies. with w_j ~ N(0, I/l^2) and
b_j ~ Uniform(0, 2*pi).

## Examples

``` r
set.seed(1L)
X_train <- matrix(rnorm(30 * 4), 30, 4)
Y_train <- matrix(rnorm(30 * 6), 30, 6)
rff <- train_rff_surrogate(X_train, Y_train, n_features = 100L)
rff$train_mse
#> [1] 0.01009735
X_new <- matrix(rnorm(5 * 4), 5, 4)
pred  <- predict_rff_surrogate(rff, X_new)
dim(pred$mean)
#> [1] 5 6
```
