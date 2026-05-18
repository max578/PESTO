# PESTO <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**P**arameter **E**stimation, **S**urrogates, and **T**ooling for **O**ptimisation

PESTO is a high-performance R package for model-independent parameter estimation and uncertainty quantification, built on modernised [PEST++](https://github.com/usgs/pestpp) algorithms with novel methodological extensions.

## Key Features

| Feature | Description |
|---------|-------------|
| **Fast ensemble solvers** | C++ (RcppEigen) implementation of IES and GLM update equations |
| **Adaptive SVD** | Automatic backend selection: randomised SVD, LAPACK, or Eigen BDCSVD |
| **Surrogate-accelerated IES** | Gaussian Process surrogates reduce model evaluations by 50--90% |
| **Adaptive ensemble sizing** | ESS-based diagnostics prevent over/under-sampling |
| **PEST++ integration** | Read/write .pst control files, run PEST++ executables from R |
| **Publication-ready plots** | Convergence, ensemble distributions, identifiability, surrogate diagnostics |

## Installation

From GitHub (development version):

```r
# install.packages("pak")
pak::pak("AAGI-AUS/PESTO")
```

Once the AAGI r-universe registry is live, pre-built binaries will
be available:

```r
install.packages("PESTO", repos = c(
  "https://aagi-aus.r-universe.dev",
  "https://cloud.r-project.org"
))
```

CRAN submission is in preparation.

### System requirements

- R >= 4.1.0
- C++17 compiler
- LAPACK/BLAS (provided by R)

## Quick Start

```r
library(PESTO)

# Create a parameter estimation scenario
pars <- data.table::data.table(
  parnme    = paste0("k", 1:10),
  partrans  = "log",
  parchglim = "factor",
  parval1   = runif(10, 0.1, 10),
  parlbnd   = 0.001,
  parubnd   = 1000,
  pargp     = "hydraulic"
)

obs <- data.table::data.table(
  obsnme = paste0("h", 1:20),
  obsval = rnorm(20, 5, 1),
  weight = 1.0,
  obgnme = "heads"
)

pst <- create_pest_scenario(pars, obs, "python model.py")
write_pst(pst, "model.pst")

# Core computational kernels (no PEST++ binary needed)
set.seed(42)
npar <- 100; nobs <- 200; nreal <- 50
upgrade <- ensemble_solution(
  par_diff  = matrix(rnorm(npar * nreal), npar, nreal),
  obs_diff  = matrix(rnorm(nobs * nreal), nobs, nreal),
  obs_resid = matrix(rnorm(nobs * nreal), nobs, nreal),
  par_resid = matrix(rnorm(npar * nreal), npar, nreal),
  weights   = abs(rnorm(nobs)) + 0.1,
  parcov_inv = abs(rnorm(npar)) + 0.1,
  Am        = matrix(rnorm(npar * (nreal - 1)), npar, nreal - 1),
  cur_lam   = 1.0
)
```

## Performance

Benchmarks on Apple Silicon (M-series) comparing PESTO's C++ kernel against a pure R implementation:

| Scale (params x obs x reals) | PESTO (ms) | R native (ms) | Speedup |
|-------------------------------|-----------|----------------|---------|
| 50 x 100 x 30                | 0.17      | 0.17           | 1.0x    |
| 500 x 1000 x 50              | 1.45      | 2.59           | 1.8x    |
| 2000 x 5000 x 100            | 22.3      | 64.0           | 2.9x    |

Randomised SVD: **35x** faster than full LAPACK at 2000x1000 for rank-20 approximation.

## Documentation

- `vignette("getting-started", package = "PESTO")` -- Introduction and basic usage
- `vignette("surrogate-ies", package = "PESTO")` -- Surrogate-accelerated IES tutorial

## Citation

```r
citation("PESTO")
```

> Moldovan, M. (2026). PESTO: Parameter Estimation, Surrogates, and Tooling for Optimisation. R package version 0.1.0. https://github.com/AAGI-AUS/PESTO

## Related Projects

- [PEST++](https://github.com/usgs/pestpp) -- The original USGS parameter estimation suite
- [pyEMU](https://github.com/pypest/pyemu) -- Python interface for PEST/PEST++

## License

GPL-3
