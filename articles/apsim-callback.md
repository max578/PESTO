# In-Process IES via R Callback — apsimx Adapter

## Why a callback driver?

[`pesto_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies.md)
shells out to the `pestpp-ies` binary, which writes and reads `.pst` /
ensemble / observation files between every realisation. For an
in-process R forward model (an `apsimx` wrapper, a Python bridge, or a
fast synthetic test problem), that file round-trip is pure overhead.
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
keeps the same C++ ensemble kernel
([`ensemble_solution()`](https://AAGI-AUS.github.io/PESTO/reference/ensemble_solution.md),
Chen & Oliver 2013) but drives the outer loop in R, calling the forward
model in-process.

This document shows:

1.  The driver on a synthetic linear-Gaussian problem (recovers truth).
2.  The shape of the
    [`apsim_callback()`](https://AAGI-AUS.github.io/PESTO/reference/apsim_callback.md)
    adapter for `apsimx` (run block is disabled — needs APSIM
    installed).

## A synthetic recovery example

``` r

library(PESTO)

npar  <- 4L
nobs  <- 8L
nreal <- 120L
sigma <- 0.05

G          <- matrix(rnorm(nobs * npar), nobs, npar)
theta_true <- c(1.0, -0.5, 2.0, 0.25)
y          <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = sigma)
names(y)   <- paste0("o", seq_len(nobs))

# Forward model: any function that takes an (nreal x npar) matrix and
# returns an (nreal x nobs) matrix.
forward <- function(theta) theta %*% t(G)

prior <- matrix(rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, paste0("p", seq_len(npar))))
```

Run six IES iterations with a single lambda:

``` r

fit <- pesto_ies_callback(
  forward_model  = forward,
  prior_ensemble = prior,
  obs            = y,
  obs_sd         = sigma,
  noptmax        = 6L,
  lambda         = 1.0,
  verbose        = FALSE
)
```

Diagnostics:

``` r

prior_rmse <- sqrt(mean((colMeans(prior) - theta_true)^2))
post_mean  <- colMeans(as.matrix(fit$par_ensemble[, -1L]))
post_rmse  <- sqrt(mean((post_mean - theta_true)^2))

cat(sprintf("prior RMSE: %.4f\n", prior_rmse))
#> prior RMSE: 1.1300
cat(sprintf("posterior RMSE: %.4f\n", post_rmse))
#> posterior RMSE: 0.0232

mean_phi <- vapply(fit$iterations, `[[`, numeric(1L), "mean_phi")
print(mean_phi)
#> [1] 17452.123014     1.534500     1.533925     1.533362     1.532812
#> [6]     1.532273
```

Posterior RMSE should be a small fraction of prior RMSE, and mean phi
should decrease across iterations.

The full per-realisation phi trajectory:

``` r

plot_phi(fit)
#> Warning in melt.data.table(phi_dt, id.vars = iter_col, measure.vars = phi_cols,
#> : 'measure.vars' [realisation, phi] are not all of the same type. By order of
#> hierarchy, the molten data value column will be of type 'double'. All measure
#> variables not of type 'double' will be coerced too. Check DETAILS in
#> ?melt.data.table for more on coercion.
```

![Per-realisation phi decreasing across IES
iterations](apsim-callback_files/figure-html/unnamed-chunk-5-1.png)

## Driving APSIM through `apsim_callback()`

The adapter wraps `apsimx` so each realisation gets its own working copy
of an APSIM template, parameter edits per `param_map`, a run, and an
extraction step. The signature returned matches what
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
expects.

``` r

library(apsimx)   # >= 2.7.0 from CRAN
# Requires a working APSIM Next Gen installation.

forward <- apsim_callback(
  template  = "wheat_wagga.apsimx",
  param_map = list(
    RUE = "Wheat.Leaf.Photosynthesis.RUE.FixedValue",
    CN2 = "Soil.SoilWater.CN2Bare"
  ),
  output_extractor = function(sim) {
    # sim is the data.frame of report-table variables.
    as.numeric(sim$Wheat.Grain.Total.Wt)
  }
)

prior <- cbind(
  RUE = runif(40, 1.0, 2.0),
  CN2 = runif(40, 60,  90)
)

fit <- pesto_ies_callback(
  forward_model  = forward,
  prior_ensemble = prior,
  obs            = c(y_2018 = 4500, y_2019 = 5200, y_2020 = 4900),
  obs_sd         = 250,
  noptmax        = 4L
)

summary(fit$par_ensemble)
```

### Failure handling

Per-realisation APSIM crashes (corrupt config, solver divergence,
missing report variable) are caught by the adapter and emerge as `NA`
rows.
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)’s
`on_failure = "na"` (default) carries those realisations forward
unchanged so the ensemble survives partial failures;
`on_failure = "stop"` aborts as soon as one row is missing.

### Concurrency

Phase-1 D4 runs realisations serially. Parallel evaluation (`future`,
`mirai`, or `apsimx`’s own job-server modes) is a planned follow-up once
`apsimx`’s thread-safety under ensemble load has been characterised. For
pure-R synthetic models the serial loop is already fast enough that the
overhead is dominated by
[`ensemble_solution()`](https://AAGI-AUS.github.io/PESTO/reference/ensemble_solution.md)
itself.

## When to prefer `pesto_ies()` instead

Use the classic `.pst` path when:

- Your forward model lives behind a non-R executable that PEST++ can
  drive directly via its template/instruction-file mechanism.
- You need full bit-for-bit compatibility with `pestpp-ies` behaviour
  (e.g. for cross-validating against an existing PEST++ workflow).
- You depend on `pestpp-ies`-specific features not yet exposed by the R
  driver (full lambda line-search, regularisation modes, Tikhonov
  priors).

Otherwise the callback driver is faster to iterate, easier to debug, and
avoids per-realisation file I/O entirely.

## Reproducibility

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] PESTO_0.3.3
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3        cli_3.6.6          knitr_1.51         rlang_1.2.0       
#>  [5] xfun_0.57          S7_0.2.2           textshaping_1.0.5  jsonlite_2.0.0    
#>  [9] data.table_1.18.4  labeling_0.4.3     glue_1.8.1         htmltools_0.5.9   
#> [13] ragg_1.5.2         sass_0.4.10        scales_1.4.0       rmarkdown_2.31    
#> [17] grid_4.6.0         evaluate_1.0.5     jquerylib_0.1.4    fastmap_1.2.0     
#> [21] yaml_2.3.12        lifecycle_1.0.5    compiler_4.6.0     RColorBrewer_1.1-3
#> [25] fs_2.1.0           Rcpp_1.1.1-1.1     farver_2.1.2       systemfonts_1.3.2 
#> [29] digest_0.6.39      R6_2.6.1           bslib_0.11.0       withr_3.0.2       
#> [33] tools_4.6.0        gtable_0.3.6       pkgdown_2.2.0      ggplot2_4.0.3     
#> [37] cachem_1.1.0       desc_1.4.3
```
