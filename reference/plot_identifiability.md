# Plot Parameter Identifiability

Creates a ranked lollipop plot of parameter identifiability based on the
singular value decomposition of the Jacobian matrix. Accepts either a
numeric Jacobian matrix in memory or a path to a `.jco` (PEST binary)
file. High-dimensional problems are kept legible by showing only the
most identifiable parameters (see `top_n`).

## Usage

``` r
plot_identifiability(
  jacobian = NULL,
  jco_file = NULL,
  pst = NULL,
  n_sv = NULL,
  top_n = NULL,
  title = "Parameter Identifiability"
)
```

## Arguments

- jacobian:

  Numeric matrix (n_obs x n_par). The Jacobian / sensitivity matrix.
  Column names, if present, are used as parameter labels; otherwise
  `p1`, `p2`, ... are generated. Either `jacobian` or `jco_file` must be
  supplied.

- jco_file:

  Character. Path to a `.jco` (Jacobian) binary file. Mutually exclusive
  with `jacobian`.

- pst:

  A `pesto_pst` object for parameter names. Optional; only used when
  reading from a `.jco` file and column names are absent.

- n_sv:

  Integer. Number of singular values to retain.

- top_n:

  Integer. Maximum number of parameters to display, ranked by
  identifiability. Defaults to all parameters when there are at most 40
  and to the 40 most identifiable otherwise, so high-dimensional
  problems stay legible. Set explicitly to override; a subtitle records
  how many of the total are shown.

- title:

  Character. Plot title.

## Value

A ggplot2 object.

## Examples

``` r
J <- matrix(rnorm(30 * 8), nrow = 30, ncol = 8)
J[, 7] <- 0.5 * J[, 1] + 0.5 * J[, 2]
J[, 8] <- 1e-6 * rnorm(30)
colnames(J) <- paste0("k", 1:8)
p <- plot_identifiability(jacobian = J)
inherits(p, "ggplot")
#> [1] TRUE
```
