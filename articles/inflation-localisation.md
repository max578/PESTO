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

A third pathology used to sit alongside these two and was by far the
largest: assimilating the *same*, unperturbed observation vector into
every realisation. That gives every member an identical data pull, so
the between-member spread stops being a posterior spread at all – and,
unlike cause 1, it does not shrink as the ensemble grows. PESTO’s
smoother now perturbs the observations per realisation under an ES-MDA
schedule (see *Posterior spread and observation perturbation* in
[`?pesto_ies_callback`](https://max578.github.io/PESTO/html/pesto_ies_callback.md)),
so what remains on this problem is the classical,
ensemble-size-dependent collapse that inflation and localisation were
designed for. The sweep below is the evidence for that: the sd-ratio now
climbs toward 1 with `nreal` instead of sitting flat.

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
    # Pin the ES-MDA observation-noise stream so every configuration below
    # is compared on identical perturbed data.
    seed           = 42L,
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

## The bare smoother under-disperses

``` r

bare <- run_ies()
round(bare$sd_ratio, 3)
#> [1] 0.916
```

The mean posterior standard deviation falls short of the analytic value
– with only 24 realisations for 6 parameters the ensemble is
over-confident. The spread-ESS ratio quantifies the same collapse from
the eigenspectrum of the parameter anomaly covariance (1 means variance
is spread isotropically across all directions; small values mean it has
collapsed onto a few):

``` r

round(bare$ess_ratio, 3)
#> [1] 0.339
```

### Is this finite-ensemble collapse, or something else?

The question is worth asking of any under-dispersed ensemble, and it has
a cheap answer. If the shortfall is finite-ensemble collapse, the
sd-ratio climbs toward 1 as the ensemble grows – a larger sample gives
the update a better estimate of the true covariance structure. If it is
flat in `nreal`, the cause is structural and no amount of inflation will
reach it. The sweep below runs the same problem at increasing `nreal`
(kept modest here for vignette build time; the pattern continues out to
8,000 realisations):

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
#> 1    24    0.897
#> 2   100    1.056
#> 3   500    1.001
```

The sd-ratio rises toward 1 across a 20-fold change in ensemble size:
the residual shortfall at `nreal = 24` is the ensemble-size-dependent
collapse, which is what inflation and localisation address. Run this
sweep on your own problem before reaching for either – a flat trend
would say the cause lies elsewhere and inflation would only be papering
over it. (This sweep is also the diagnostic that exposed the
unperturbed-observation defect: before PESTO perturbed the observations
the same sweep sat flat at about 0.21 from 24 realisations to 8,000.)

## Inflation re-expands the spread

[`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)
offers four methods. The workhorse is relaxation to prior spread
(`"rtps"`, Whitaker & Hamill 2012): each parameter’s posterior anomalies
are rescaled toward the pre-update spread, so the directions that
collapsed hardest are re-inflated most. The `"adaptive"` method instead
targets a global spread-retention floor.

``` r

rtps     <- run_ies(inflation = pesto_inflation("rtps", alpha = 0.2))
adaptive <- run_ies(inflation = pesto_inflation("adaptive",
                                                retention_floor = 0.6))

data.frame(
  method        = c("none", "rtps", "adaptive"),
  sd_ratio      = round(c(bare$sd_ratio, rtps$sd_ratio, adaptive$sd_ratio), 3),
  spread_ess    = round(c(bare$ess_ratio, rtps$ess_ratio, adaptive$ess_ratio), 3)
)
#>     method sd_ratio spread_ess
#> 1     none    0.916      0.339
#> 2     rtps    1.060      0.344
#> 3 adaptive    0.947      0.319
```

Both methods lift the retained posterior spread toward the analytic
value. The caveat worth stating plainly is the other one: inflation is a
*tunable* remedy, not a self-correcting one. The strengths used here
(`alpha = 0.2`, `retention_floor = 0.6`) were chosen against this
problem’s known posterior; a stronger setting overshoots and returns an
ensemble that is too wide, which is no more honest than one that is too
narrow. Measure the sd-ratio or the interval coverage on your own
problem and tune to it – do not carry these numbers across.

The spread-ESS diagnostic
([`ensemble_spread_ess()`](https://max578.github.io/PESTO/html/ensemble_spread_ess.md))
is recorded on **every** iteration regardless of method, so the collapse
trajectory is always available in the result:

``` r

ess_trace <- function(inflation = NULL) {
  fit <- pesto_ies_callback(
    forward, prior, setNames(y, paste0("o", seq_len(nobs))),
    obs_sd = obs_sd, noptmax = 12L, inflation = inflation,
    seed = 42L, verbose = FALSE
  )
  vapply(fit$iterations, function(d) d$spread_ess_ratio, numeric(1L))
}
tr_bare <- ess_trace()
tr_rtps <- ess_trace(pesto_inflation("rtps", alpha = 0.2))

plot(tr_bare, type = "b", pch = 19, col = "grey40", ylim = c(0, 1),
     xlab = "iteration", ylab = "spread-ESS ratio",
     main = "Dispersion trajectory, bare vs RTPS inflation")
lines(tr_rtps, type = "b", pch = 17, col = "#1B7837")
legend("topright", c("bare", "RTPS"), pch = c(19, 17),
       col = c("grey40", "#1B7837"), bty = "n")
```

![Two lines plotted against IES iteration number (x axis, 1 to 12) with
the spread-ESS ratio on the y axis from 0 to 1. Both start near 0.5 to
0.6 and decline gradually to between 0.3 and 0.35 by iteration 12. The
inflated line sits consistently above the bare line by roughly
0.04.](inflation-localisation_files/figure-html/ess-trace-1.png)

Spread-ESS ratio by IES iteration, without and with RTPS inflation. Both
decline gently as the update concentrates the ensemble on the informed
directions; the inflated run sits above the bare one at every iteration,
which is what inflation buys here – a slower drain, not a flat line.

Read the gap between the two lines, not the shape of either. The gentle
decline is the update doing its job – concentrating the ensemble on the
directions the data inform – and the inflated run holding above the bare
one at every iteration is the effect being demonstrated. The trace is a
trajectory diagnostic; the sd-ratio numbers above are the calibration
check, comparing the *final* ensemble against the analytic posterior.

## Localisation suppresses spurious correlations

For parameter-estimation problems whose parameters carry no spatial
coordinate, the recommended localiser is the correlation-based automatic
method (Luo & Bhakta 2020). It needs no metric: it estimates a noise
floor from the ensemble itself and damps sample correlations that fall
below it.

``` r

loc <- run_ies(localisation = pesto_localisation("correlation",
                                                 taper = "hard"))
round(loc$sd_ratio, 3)
#> [1] 0.954
```

The two countermeasures compose – inflation restores variance magnitude,
localisation removes the spurious updates that drain it:

``` r

both <- run_ies(
  inflation    = pesto_inflation("rtps", alpha = 0.1),
  localisation = pesto_localisation("correlation", taper = "hard")
)
round(both$sd_ratio, 3)
#> [1] 1.045
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
