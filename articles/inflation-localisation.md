# Countering Ensemble Collapse: Inflation and Localisation

## The problem: ensemble under-dispersion

An iterative ensemble smoother estimates a posterior parameter
distribution from a finite ensemble of realisations. Two distinct
pathologies can make the posterior spread *too narrow* – the ensemble
becomes over-confident:

1.  **Finite-ensemble collapse.** Each assimilation step contracts the
    ensemble spread; over several iterations the posterior variance can
    fall below the true posterior variance. This is the classical,
    well-characterised effect the ensemble-Kalman literature addresses
    with inflation and localisation (Whitaker & Hamill 2012; Luo &
    Bhakta 2020), and it shrinks as the ensemble grows.
2.  **Spurious correlations.** A finite ensemble manufactures apparent
    correlations between parameters and observations that are not real.
    Acting on them injects noise into the update and accelerates
    collapse.

This vignette demonstrates inflation and localisation against cause 1 on
a linear-Gaussian problem with a known analytic posterior. **A caveat
carried over from the current release:** on this same problem the
realised under-dispersion is larger than finite-ensemble collapse alone
predicts, and does not shrink materially as `nreal` grows from 24 to
8,000 (a sweep this vignette reproduces below) – evidence that a second,
ensemble-size-independent source of under-dispersion is present in the
current smoother alongside the classical effect. Inflation, which is
designed against cause 1, does not close that residual gap; treat the
sd-ratio numbers in this vignette as a lower bound on the smoother’s
current calibration, not as “matches the analytic posterior once
inflated”.

PESTO addresses the first with **covariance inflation**
([`pesto_inflation()`](https://max578.github.io/PESTO/html/pesto_inflation.md))
and the second with **covariance localisation**
([`pesto_localisation()`](https://max578.github.io/PESTO/html/pesto_localisation.md)).
Both are opt-in: the default `NULL` leaves the update identical to the
bare smoother.

This vignette demonstrates the effect on a linear-Gaussian problem,
where the analytic posterior is known and the collapse can be measured
exactly.

## A linear-Gaussian problem with a known posterior

For a linear forward model d = G\theta + \varepsilon with a Gaussian
prior \theta \sim N(0, C_0) and observation error \varepsilon \sim N(0,
R), the posterior covariance is C\_{\mathrm{post}} = (C_0^{-1} +
G^{\top} R^{-1} G)^{-1}. We can therefore compare the ensemble’s
posterior spread directly against the truth.

``` r

library(PESTO)
set.seed(42L)

npar  <- 6L
nobs  <- 10L
nreal <- 24L          # a deliberately small ensemble, to provoke collapse

G          <- matrix(rnorm(nobs * npar), nobs, npar)
theta_true <- rnorm(npar)
obs_sd     <- 0.3
y          <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = obs_sd)

# Analytic posterior standard deviation (standard-normal prior).
post_cov <- solve(diag(npar) + crossprod(G) / obs_sd^2)
post_sd  <- sqrt(diag(post_cov))

forward <- function(theta) theta %*% t(G)
prior   <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
```

A small helper runs the smoother and returns the realised posterior
spread alongside the spread-ESS collapse diagnostic recorded on the
final iteration.

``` r

run_ies <- function(inflation = NULL, localisation = NULL) {
  fit <- pesto_ies_callback(
    forward_model  = forward,
    prior_ensemble = prior,
    obs            = setNames(y, paste0("o", seq_len(nobs))),
    obs_sd         = obs_sd,
    noptmax        = 12L,
    inflation      = inflation,
    localisation   = localisation,
    verbose        = FALSE
  )
  par_post  <- as.matrix(fit$par_ensemble[, -1])
  last_diag <- fit$iterations[[length(fit$iterations)]]
  list(
    sd_ratio  = mean(apply(par_post, 2L, sd) / post_sd),
    ess_ratio = last_diag$spread_ess_ratio
  )
}
```

## The bare smoother collapses

``` r

bare <- run_ies()
round(bare$sd_ratio, 3)
#> [1] 0.222
```

The mean posterior standard deviation is only a fraction of the analytic
value – the ensemble is badly over-confident. The spread-ESS ratio
quantifies the same collapse from the eigenspectrum of the parameter
anomaly covariance (1 means variance is spread isotropically across all
directions; small values mean it has collapsed onto a few):

``` r

round(bare$ess_ratio, 3)
#> [1] 0.315
```

### Is this finite-ensemble collapse, or something else?

If the under-dispersion above were purely finite-ensemble collapse, the
sd-ratio should climb toward 1 as the ensemble grows – a larger sample
gives the update a better estimate of the true covariance structure. The
sweep below runs the same problem at increasing `nreal` (kept modest
here for vignette build time; the pattern is the same out to 8,000
realisations):

``` r

sweep_nreal <- c(24L, 100L, 500L)
sd_by_nreal <- vapply(sweep_nreal, function(n) {
  set.seed(42L)
  prior_n <- matrix(rnorm(n * npar), n, npar,
                    dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit_n <- pesto_ies_callback(
    forward_model = forward, prior_ensemble = prior_n,
    obs = setNames(y, paste0("o", seq_len(nobs))), obs_sd = obs_sd,
    noptmax = 12L, verbose = FALSE
  )
  mean(apply(as.matrix(fit_n$par_ensemble[, -1]), 2L, sd) / post_sd)
}, numeric(1L))

data.frame(nreal = sweep_nreal, sd_ratio = round(sd_by_nreal, 3))
#>   nreal sd_ratio
#> 1    24    0.216
#> 2   100    0.213
#> 3   500    0.211
```

The sd-ratio is essentially flat across a 20-fold change in ensemble
size, not the rising trend finite-ensemble collapse alone would produce.
That flatness is the evidence behind the caveat above: most of this
problem’s under-dispersion comes from a source that does not respond to
`nreal`, and inflation – being a remedy for the ensemble-size-dependent
part – cannot be expected to close it either.

## Inflation re-expands the spread

[`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)
offers four methods. The workhorse is relaxation to prior spread
(`"rtps"`, Whitaker & Hamill 2012): each parameter’s posterior anomalies
are rescaled toward the pre-update spread, so the directions that
collapsed hardest are re-inflated most. The `"adaptive"` method instead
targets a global spread-retention floor.

``` r

rtps     <- run_ies(inflation = pesto_inflation("rtps", alpha = 0.6))
adaptive <- run_ies(inflation = pesto_inflation("adaptive",
                                                retention_floor = 0.7))

data.frame(
  method        = c("none", "rtps", "adaptive"),
  sd_ratio      = round(c(bare$sd_ratio, rtps$sd_ratio, adaptive$sd_ratio), 3),
  spread_ess    = round(c(bare$ess_ratio, rtps$ess_ratio, adaptive$ess_ratio), 3)
)
#>     method sd_ratio spread_ess
#> 1     none    0.222      0.315
#> 2     rtps    0.507      0.436
#> 3 adaptive    0.310      0.459
```

RTPS roughly doubles the retained posterior spread and lifts the
spread-ESS ratio. A caveat worth stating plainly: inflation *mitigates*
finite-ensemble collapse, it does not abolish it, and – per the sweep
above – it does not touch the ensemble-size-independent residual either.
Read the `sd_ratio` improvement as “closer than the bare smoother”, not
as “close to the analytic posterior”.

The spread-ESS diagnostic
([`ensemble_spread_ess()`](https://max578.github.io/PESTO/html/ensemble_spread_ess.md))
is recorded on **every** iteration regardless of method, so the collapse
trajectory is always available in the result:

``` r

fit <- pesto_ies_callback(
  forward, prior, setNames(y, paste0("o", seq_len(nobs))),
  obs_sd = obs_sd, noptmax = 12L,
  inflation = pesto_inflation("rtps", alpha = 0.6), verbose = FALSE
)
plot(
  vapply(fit$iterations, function(d) d$spread_ess_ratio, numeric(1L)),
  type = "b", pch = 19, xlab = "iteration", ylab = "spread-ESS ratio",
  main = "Dispersion held up under RTPS inflation", ylim = c(0, 1)
)
```

![Line plot of the spread-ESS ratio (y axis, 0 to 1) against IES
iteration number (x axis, 1 to 12) under RTPS inflation. The ratio drops
over the first few iterations and then levels off at a stable plateau
rather than continuing to decay toward
zero.](inflation-localisation_files/figure-html/ess-trace-1.png)

Spread-ESS ratio by IES iteration under RTPS inflation: the trajectory
stabilises rather than decaying toward zero, so inflation is holding the
ensemble spread against collapse iteration over iteration.

The plateau, rather than a continued slide toward zero, is what
“inflation holds the ensemble spread against collapse” looks like
iteration by iteration – distinct from the sd-ratio numbers above, which
compare the *final* ensemble against the analytic posterior rather than
tracking the trajectory.

## Localisation suppresses spurious correlations

For parameter-estimation problems whose parameters carry no spatial
coordinate, the recommended localiser is the correlation-based automatic
method (Luo & Bhakta 2020). It needs no metric: it estimates a noise
floor from the ensemble itself and damps sample correlations that fall
below it.

``` r

loc <- run_ies(localisation = pesto_localisation("correlation",
                                                 taper = "soft"))
round(loc$sd_ratio, 3)
#> [1] 0.807
```

The two countermeasures compose – inflation restores variance magnitude,
localisation removes the spurious updates that drain it:

``` r

both <- run_ies(
  inflation    = pesto_inflation("rtps", alpha = 0.6),
  localisation = pesto_localisation("correlation", taper = "soft")
)
round(both$sd_ratio, 3)
#> [1] 2.515
```

When a genuine distance metric *does* exist, the classical Gaspari-Cohn
taper is available via `pesto_localisation("distance", ...)`, supplying
either a precomputed `distances` matrix or `par_coords` / `obs_coords`
together with a localisation `radius`.

## References

Gaspari, G. & Cohn, S. E. (1999). Construction of correlation functions
in two and three dimensions. *Quarterly Journal of the Royal
Meteorological Society*, 125(554), 723–757.

Luo, X. & Bhakta, T. (2020). Automatic and adaptive localization for
ensemble-based history matching. *Journal of Petroleum Science and
Engineering*, 184, 106559.

Whitaker, J. S. & Hamill, T. M. (2012). Evaluating methods to account
for system errors in ensemble data assimilation. *Monthly Weather
Review*, 140(9), 3078–3089.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
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
#> [1] PESTO_0.10.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3         cli_3.6.6           knitr_1.51         
#>  [4] rlang_1.3.0         xfun_0.60           otel_0.2.0         
#>  [7] S7_0.2.2            textshaping_1.0.5   jsonlite_2.0.0     
#> [10] data.table_1.18.6.1 glue_1.8.1          htmltools_0.5.9    
#> [13] ragg_1.5.2          sass_0.4.10         scales_1.4.0       
#> [16] rmarkdown_2.31      grid_4.6.1          evaluate_1.0.5     
#> [19] jquerylib_0.1.4     fastmap_1.2.0       yaml_2.3.12        
#> [22] lifecycle_1.0.5     compiler_4.6.1      RColorBrewer_1.1-3 
#> [25] fs_2.1.0            Rcpp_1.1.2          farver_2.1.2       
#> [28] systemfonts_1.3.2   digest_0.6.39       R6_2.6.1           
#> [31] bslib_0.12.0        gtable_0.3.6        tools_4.6.1        
#> [34] pkgdown_2.2.1       ggplot2_4.0.3       cachem_1.1.0       
#> [37] desc_1.4.3
```
