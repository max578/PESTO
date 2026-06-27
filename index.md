# PESTO ![](reference/figures/logo.png)

**P**arameter **E**stimation, **S**urrogates, and **T**ooling for
**O**ptimisation

PESTO is a high-performance R package for model-independent parameter
estimation and uncertainty quantification. It brings the algorithms of
**PEST** (*Parameter ESTimation*; Doherty 2015) and its C++ successor
[**PEST++**](https://github.com/usgs/pestpp) (White et al. 2020) –
notably the iterative ensemble smoother (IES) – natively into R, and is
the first R-native implementation of that algorithm family. It adds a
typed forward-model contract, an in-process simulator callback,
multi-fidelity acceleration, covariance inflation / localisation, and
surrogate methods.

## Key Features

| Feature | Description |
|----|----|
| **Fast ensemble solvers** | C++ (RcppEigen) implementation of IES and GLM update equations |
| **In-process callback** | Couple any R forward model directly – no file-exchange overhead |
| **Adaptive SVD** | Automatic backend selection: randomised SVD, LAPACK, or Eigen BDCSVD |
| **Surrogate-accelerated IES** | Gaussian Process surrogates reduce model evaluations by 50–90% |
| **Inflation & localisation** | Counter finite-ensemble under-dispersion and spurious correlations |
| **Adaptive ensemble sizing** | ESS-based diagnostics prevent over/under-sampling |
| **PEST++ integration** | Read/write .pst control files, run PEST++ executables from R |
| **Publication-ready plots** | Convergence, ensemble distributions, identifiability, surrogate diagnostics |

## Installation

The canonical home for PESTO is
[`max578/PESTO`](https://github.com/max578/PESTO), published through the
author’s personal CRAN-track channel and the `max578.r-universe.dev`
registry. The `AAGI-AUS/PESTO` remote is retained as a read-only mirror.

``` r

# Development version from GitHub
# install.packages("pak")
pak::pak("max578/PESTO")

# Pre-built binaries from r-universe
install.packages("PESTO", repos = c(
  "https://max578.r-universe.dev",
  "https://cloud.r-project.org"
))
```

CRAN submission is in preparation.

## Ensemble Bayesian inference (IES)

PESTO performs approximate Bayesian inference for parameter estimation.
You specify a **prior as a parameter ensemble** – a matrix of draws
encoding your prior belief about each parameter (rows are realisations,
columns are parameters) – together with the observations and their error
standard deviation. PESTO conditions the ensemble on the observations
through the IES update and returns an approximate **posterior ensemble**
with uncertainty diagnostics. Any R function mapping a parameter matrix
to an observation matrix can be the forward model.

``` r

library(PESTO)
set.seed(1)

# Forward model: y = G %*% theta (any R function theta-matrix -> obs-matrix works)
G <- matrix(rnorm(6 * 3), 6, 3)
forward <- function(theta) theta %*% t(G)
truth <- c(a = 1, b = -0.5, c = 2)
obs <- stats::setNames(as.numeric(G %*% truth) + rnorm(6, sd = 0.05),
                       paste0("o", 1:6))

# Prior: an ensemble of parameter draws (here a broad Gaussian over 3 parameters)
prior <- matrix(rnorm(80 * 3), 80, 3, dimnames = list(NULL, names(truth)))

fit <- pesto_ies_callback(forward, prior, obs, obs_sd = 0.05, noptmax = 6,
                          verbose = FALSE)
colMeans(as.matrix(fit$par_ensemble[, -1]))   # posterior mean -> (1, -0.5, 2)
#>          a          b          c 
#>  1.0136949 -0.4936717  1.9847292
```

Optional convergence-based early stopping (`phi_tol`), covariance
inflation
([`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)),
and localisation
([`pesto_localisation()`](https://max578.github.io/PESTO/reference/pesto_localisation.md))
are documented in
[`?pesto_ies_callback`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
and the *Getting started* vignette.

## Low-level kernels

The C++ update kernels are exported for advanced use (no PEST++ binary
needed):

``` r

set.seed(42)
npar <- 100; nobs <- 200; nreal <- 50
upgrade <- ensemble_solution(
  par_diff   = matrix(rnorm(npar * nreal), npar, nreal),
  obs_diff   = matrix(rnorm(nobs * nreal), nobs, nreal),
  obs_resid  = matrix(rnorm(nobs * nreal), nobs, nreal),
  par_resid  = matrix(rnorm(npar * nreal), npar, nreal),
  weights    = abs(rnorm(nobs)) + 0.1,
  parcov_inv = abs(rnorm(npar)) + 0.1,
  Am         = matrix(rnorm(npar * (nreal - 1)), npar, nreal - 1),
  cur_lam    = 1.0
)
dim(upgrade)
#> [1]  50 100
```

The classical `.pst`-file path is also supported:

``` r

pars <- data.table::data.table(
  parnme = paste0("k", 1:10), partrans = "log", parchglim = "factor",
  parval1 = runif(10, 0.1, 10), parlbnd = 0.001, parubnd = 1000,
  pargp = "hydraulic"
)
obs <- data.table::data.table(
  obsnme = paste0("h", 1:20), obsval = rnorm(20, 5, 1),
  weight = 1.0, obgnme = "heads"
)
pst <- create_pest_scenario(pars, obs, "python model.py")
write_pst(pst, file.path(tempdir(), "model.pst"))
```

## Dependencies

- R \>= 4.1.0
- C++17 compiler (clang on macOS, g++ on Linux, Rtools on Windows)
- LAPACK/BLAS (bundled with R)

R-package dependencies are declared in `DESCRIPTION` (Imports: `Rcpp`,
`data.table`, `ggplot2`, `S7`, `yaml`, `digest`).

## Documentation

- [`vignette("getting-started", package = "PESTO")`](https://max578.github.io/PESTO/articles/getting-started.md)
  – Introduction and basic usage
- [`vignette("inflation-localisation", package = "PESTO")`](https://max578.github.io/PESTO/articles/inflation-localisation.md)
  – Finite-ensemble countermeasures
- [`vignette("surrogate-ies", package = "PESTO")`](https://max578.github.io/PESTO/articles/surrogate-ies.md)
  – Surrogate-accelerated IES
- [`vignette("pestpp-comparison-and-simulation", package = "PESTO")`](https://max578.github.io/PESTO/articles/pestpp-comparison-and-simulation.md)
  – Benchmark vs PEST/PEST++
- [`vignette("ensemble-manifest", package = "PESTO")`](https://max578.github.io/PESTO/articles/ensemble-manifest.md)
  – The reproducible run manifest

## Contributing

Contributions are welcome. See
[`CONTRIBUTING.md`](https://max578.github.io/PESTO/CONTRIBUTING.md) for
the development workflow and pull-request convention; bug reports and
feature requests go through [GitHub
Issues](https://github.com/max578/PESTO/issues). All participants abide
by the [Code of
Conduct](https://max578.github.io/PESTO/CODE_OF_CONDUCT.md).

## Citation

``` r

citation("PESTO")
```

> Moldovan, M. (2026). PESTO: Parameter Estimation, Surrogates, and
> Tooling for Optimisation. R package version 0.8.0.
> <https://github.com/max578/PESTO>

## Acknowledgements

PESTO builds on the algorithmic legacy of the
[PEST++](https://github.com/usgs/pestpp) project (US Geological Survey)
and the underlying PEST framework by John Doherty. Developed at Adelaide
University; the surrogate-acceleration, adaptive-ensemble-sizing,
multi-fidelity, and convergence-aware components are original
contributions of the author.

## License

GPL (\>= 3). See
[`LICENSE.md`](https://max578.github.io/PESTO/LICENSE.md) for the full
text.
