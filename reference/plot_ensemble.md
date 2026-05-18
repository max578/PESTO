# Plot Ensemble Parameter Distributions

Visualises the prior and/or posterior parameter ensemble distributions
as violin plots or density ridges.

## Usage

``` r
plot_ensemble(
  ensemble,
  parameters = NULL,
  prior_ensemble = NULL,
  max_params = 20L,
  title = "Parameter Ensemble Distributions"
)
```

## Arguments

- ensemble:

  data.table. Parameter ensemble (rows = realisations).

- parameters:

  Character vector. Parameter names to plot. If NULL, selects up to 20
  parameters with highest variance.

- prior_ensemble:

  data.table. Optional prior ensemble for comparison.

- max_params:

  Integer. Maximum parameters to display.

- title:

  Character. Plot title.

## Value

A ggplot2 object.

## Examples

``` r
posterior <- data.table::data.table(
  real_name = sprintf("real_%02d", 1:50),
  k1 = rnorm(50, 1.0, 0.15),
  k2 = rnorm(50, 0.5, 0.05),
  k3 = rnorm(50, 2.0, 0.40),
  k4 = rnorm(50, 0.1, 0.02)
)
prior <- data.table::data.table(
  real_name = sprintf("real_%02d", 1:50),
  k1 = rnorm(50, 1.0, 0.50),
  k2 = rnorm(50, 0.5, 0.20),
  k3 = rnorm(50, 2.0, 1.00),
  k4 = rnorm(50, 0.1, 0.08)
)
p <- plot_ensemble(posterior, prior_ensemble = prior)
inherits(p, "ggplot")
#> [1] TRUE
```
