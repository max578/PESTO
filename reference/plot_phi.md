# Plot Objective Function (Phi) Convergence

Creates a publication-quality plot of objective function values across
iterations, showing mean, min, max, and individual realisation traces.

## Usage

``` r
plot_phi(
  result,
  log_scale = TRUE,
  show_reals = FALSE,
  title = "Objective Function Convergence"
)
```

## Arguments

- result:

  A `pesto_ies_result` or `pesto_glm_result` object, or a data.table
  with columns `iteration` and `phi` (at minimum).

- log_scale:

  Logical. Use log10 scale for y-axis (default TRUE).

- show_reals:

  Logical. Show individual realisation traces.

- title:

  Character. Plot title.

## Value

A ggplot2 object.

## Examples

``` r
phi_dt <- data.table::data.table(
  iteration = 0:4,
  total_runs = c(50L, 100L, 150L, 200L, 250L),
  mean = c(1200, 450, 180, 95, 72),
  min  = c(900, 320, 130, 70, 55),
  max  = c(1700, 680, 260, 140, 105),
  median = c(1180, 440, 175, 92, 70),
  std    = c(180, 80, 35, 18, 12)
)
p <- plot_phi(phi_dt, log_scale = TRUE,
              title = "Synthetic Phi Convergence")
inherits(p, "ggplot")
#> [1] TRUE
```
