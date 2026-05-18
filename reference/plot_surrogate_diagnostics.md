# Plot Surrogate Diagnostics

Visualises the surrogate-accelerated IES performance including model
savings, uncertainty distribution, and GP quality metrics.

## Usage

``` r
plot_surrogate_diagnostics(results, title = "Surrogate IES Diagnostics")
```

## Arguments

- results:

  List of surrogate update results from multiple iterations.

- title:

  Character. Plot title.

## Value

A ggplot2 object.

## Examples

``` r
iter1 <- list(n_model_runs = 12L, n_surrogate_runs = 38L,
              savings_pct = 76.0, mean_uncertainty = 0.18)
iter2 <- list(n_model_runs = 8L,  n_surrogate_runs = 42L,
              savings_pct = 84.0, mean_uncertainty = 0.11)
iter3 <- list(n_model_runs = 5L,  n_surrogate_runs = 45L,
              savings_pct = 90.0, mean_uncertainty = 0.07)
p <- plot_surrogate_diagnostics(list(iter1, iter2, iter3))
inherits(p, "ggplot")
#> [1] TRUE
```
