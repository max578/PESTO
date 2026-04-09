# PESTO 0.1.0

## Initial Release

### Core Features
* `ensemble_solution()` — High-performance C++ implementation of the IES
  ensemble update equation (Chen & Oliver, 2013) via RcppEigen.
* `ensemble_solution_mda()` — Multiple Data Assimilation (Evensen, 2018)
  update kernel.
* `compute_phi()` — Fast weighted sum-of-squares objective function.

### Hardware Acceleration
* `adaptive_svd()` — Automatic SVD backend selection (LAPACK, Eigen BDCSVD,
  or randomised SVD) based on matrix size and target rank.
* `rsvd()` — Randomised SVD (Halko-Martinsson-Tropp, 2011) for asymptotically
  faster rank-k approximations.
* `accelerate_svd()` — Direct LAPACK SVD leveraging platform-optimised BLAS
  (Apple Accelerate/AMX on macOS, MKL or OpenBLAS on Linux).
* `ensemble_solution_gpu()` — GPU-ready ensemble solution with adaptive SVD
  backend and performance diagnostics.

### Surrogate-Accelerated IES (Novel)
* `train_gp_surrogate()` — Gaussian Process surrogate model training with
  automatic hyperparameter selection (median heuristic).
* `predict_gp_surrogate()` — GP prediction with uncertainty quantification.
* `surrogate_ensemble_update()` — Surrogate-accelerated IES update with
  adaptive model/surrogate switching and control-variate bias correction.
* `adaptive_ensemble_size()` — Convergence-aware dynamic ensemble sizing
  based on ESS and coefficient of variation diagnostics.

### PEST++ Integration
* `read_pst()` / `write_pst()` — PEST control file I/O.
* `read_ensemble()` / `write_ensemble()` — Ensemble file I/O (CSV + binary).
* `pesto_ies()`, `pesto_glm()`, `pesto_sweep()`, `pesto_sensitivity()` —
  High-level wrappers for PEST++ executables.
* `create_pest_scenario()` — Programmatic scenario builder.

### Visualisation
* `plot_phi()` — Objective function convergence plotting.
* `plot_ensemble()` — Prior/posterior parameter distribution comparison.
* `plot_identifiability()` — SVD-based parameter identifiability analysis.
* `plot_surrogate_diagnostics()` — Surrogate IES performance visualisation.
