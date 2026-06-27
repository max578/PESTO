# Package index

## Core ensemble solvers

C++ kernel update equations (Chen-Oliver IES, Evensen MDA) and weighted
objective.

- [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  : Ensemble Solution Kernel (GLM form)
- [`ensemble_solution_mda()`](https://max578.github.io/PESTO/reference/ensemble_solution_mda.md)
  : Ensemble Solution Kernel (MDA / Evensen form)
- [`ensemble_solution_gpu()`](https://max578.github.io/PESTO/reference/ensemble_solution_gpu.md)
  : Ensemble Solution with Adaptive SVD Backend
- [`compute_phi()`](https://max578.github.io/PESTO/reference/compute_phi.md)
  : Compute Phi (Objective Function) for Ensemble

## SVD backends

Automatic backend selection; randomised SVD; direct LAPACK.

- [`adaptive_svd()`](https://max578.github.io/PESTO/reference/adaptive_svd.md)
  : Adaptive SVD with Automatic Backend Selection
- [`accelerate_svd()`](https://max578.github.io/PESTO/reference/accelerate_svd.md)
  : Hardware-Accelerated SVD via LAPACK
- [`rsvd()`](https://max578.github.io/PESTO/reference/rsvd.md) :
  Randomised SVD (Halko-Martinsson-Tropp Algorithm)

## High-level PEST++ wrappers

Convenience wrappers driving PEST++ executables (`pestpp-ies`,
`pestpp-glm`, …) from R.

- [`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md)
  : Run PEST++ IES (Iterative Ensemble Smoother)
- [`pesto_glm()`](https://max578.github.io/PESTO/reference/pesto_glm.md)
  : Run PEST++ GLM (Gauss-Levenberg-Marquardt)
- [`pesto_sweep()`](https://max578.github.io/PESTO/reference/pesto_sweep.md)
  : Run PEST++ SWP (Parametric Sweep)
- [`pesto_sensitivity()`](https://max578.github.io/PESTO/reference/pesto_sensitivity.md)
  : Run PEST++ SEN (Global Sensitivity Analysis)
- [`pesto_surrogate_ies()`](https://max578.github.io/PESTO/reference/pesto_surrogate_ies.md)
  : Surrogate-Accelerated IES Iteration
- [`pesto_version()`](https://max578.github.io/PESTO/reference/pesto_version.md)
  : Get PESTO package version information

## Forward-model contract and multi-fidelity

Typed forward-model object both adapter modes honour, and the
multi-fidelity (cheap / expensive) bridge with its control-variate
combiner.

- [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  : Forward-Model Contract (S7 class)

- [`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md)
  :

  Coerce an object into a `pesto_forward_model`

- [`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md)
  : Evaluate a PESTO forward model

- [`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
  : Multi-Fidelity Forward Model (S7 class)

- [`mf_control_variate()`](https://max578.github.io/PESTO/reference/mf_control_variate.md)
  : Affine control-variate bias correction across fidelities

## ODE / compartmental forward-model templates

Ready-to-use differential-equation forward models – a generic ODE
builder plus crop-growth and SEIR specialisations – that plug into the
IES driver as typed forward-model objects.

- [`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
  : Forward Model from a System of Ordinary Differential Equations
- [`crop_growth_forward_model()`](https://max578.github.io/PESTO/reference/crop_growth_forward_model.md)
  : Crop-Growth Forward Model (Logistic / Expolinear Biomass)
- [`seir_forward_model()`](https://max578.github.io/PESTO/reference/seir_forward_model.md)
  : SEIR Compartmental Forward Model

## In-process IES via R callback

Driver and APSIM adapter for forward-model callable IES, bypassing the
.pst file cycle.

- [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  : Run IES with an In-Process R Callback Forward Model
- [`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md)
  : Run a Sequential (Filter-Mode) Iterative Ensemble Smoother
- [`print(`*`<pesto_ies_result>`*`)`](https://max578.github.io/PESTO/reference/pesto_ies_result-methods.md)
  [`plot(`*`<pesto_ies_result>`*`)`](https://max578.github.io/PESTO/reference/pesto_ies_result-methods.md)
  : Print and plot methods for PESTO ensemble-smoother results
- [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  : apsimx Forward-Model Adapter for PESTO IES

## Inflation and localisation

Finite-ensemble pathology countermeasures – covariance inflation against
under-dispersion, localisation against spurious correlations, and the
spread-ESS collapse diagnostic.

- [`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md)
  : Covariance Inflation Specification for IES
- [`pesto_localisation()`](https://max578.github.io/PESTO/reference/pesto_localisation.md)
  : Covariance Localisation Specification for IES
- [`ensemble_spread_ess()`](https://max578.github.io/PESTO/reference/ensemble_spread_ess.md)
  : Spectral Spread Effective Sample Size of a Parameter Ensemble
- [`correlation_localisation()`](https://max578.github.io/PESTO/reference/correlation_localisation.md)
  : Correlation-Based Automatic Localisation Taper
- [`gaspari_cohn()`](https://max578.github.io/PESTO/reference/gaspari_cohn.md)
  : Gaspari-Cohn Localisation Taper
- [`ensemble_solution_localised()`](https://max578.github.io/PESTO/reference/ensemble_solution_localised.md)
  : Localised Ensemble Solution Kernel (explicit-gain GLM form)

## Reference IES (pure R)

Textbook Chen & Oliver (2013) reference implementation used as the
canonical comparison target by the comparison vignette.

- [`pesto_reference_ies()`](https://max578.github.io/PESTO/reference/pesto_reference_ies.md)
  : Pure-R Reference IES Update — Chen & Oliver (2013) eq. 12

## Surrogates

Gaussian Process and Random Fourier Feature surrogates;
convergence-aware ensemble sizing.

- [`train_gp_surrogate()`](https://max578.github.io/PESTO/reference/train_gp_surrogate.md)
  : Train a Gaussian Process Surrogate
- [`predict_gp_surrogate()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate.md)
  : Predict with GP Surrogate (with Uncertainty)
- [`train_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/train_gp_surrogate_tuned.md)
  : Train a GP Surrogate with Maximum-Likelihood (Anisotropic) Length
  Scales
- [`predict_gp_surrogate_tuned()`](https://max578.github.io/PESTO/reference/predict_gp_surrogate_tuned.md)
  : Predict from an MLE-Tuned GP Surrogate
- [`train_rff_surrogate()`](https://max578.github.io/PESTO/reference/train_rff_surrogate.md)
  : Train a Sparse GP Surrogate via Random Fourier Features
- [`predict_rff_surrogate()`](https://max578.github.io/PESTO/reference/predict_rff_surrogate.md)
  : Predict with RFF Sparse GP Surrogate
- [`surrogate_ensemble_update()`](https://max578.github.io/PESTO/reference/surrogate_ensemble_update.md)
  : Surrogate-Accelerated Ensemble Update
- [`adaptive_ensemble_size()`](https://max578.github.io/PESTO/reference/adaptive_ensemble_size.md)
  : Adaptive Ensemble Sizing
- [`check_surrogate_regime()`](https://max578.github.io/PESTO/reference/check_surrogate_regime.md)
  : Check whether a surrogate-IES regime is favourable

## PEST control file I/O

Read/write `.pst` files and ensemble snapshots; build scenarios
programmatically.

- [`read_pst()`](https://max578.github.io/PESTO/reference/read_pst.md) :
  Read a PEST Control File (.pst)
- [`write_pst()`](https://max578.github.io/PESTO/reference/write_pst.md)
  : Write a PEST Control File (.pst)
- [`create_pest_scenario()`](https://max578.github.io/PESTO/reference/create_pest_scenario.md)
  : Create a PEST Scenario Programmatically
- [`read_ensemble()`](https://max578.github.io/PESTO/reference/read_ensemble.md)
  : Read an Ensemble File
- [`write_ensemble()`](https://max578.github.io/PESTO/reference/write_ensemble.md)
  : Write an Ensemble File
- [`print(`*`<pesto_pst>`*`)`](https://max578.github.io/PESTO/reference/print.pesto_pst.md)
  : Print method for pesto_pst objects

## Ensemble manifest (S7 cross-package contract)

Versioned, hashed, provenance-tracked ensemble-run container consumed by
kernR and proxymix.

- [`pesto_ensemble_manifest`](https://max578.github.io/PESTO/reference/pesto_ensemble_manifest.md)
  : PESTO Ensemble Manifest (S7 class)

- [`pesto_obs_schema()`](https://max578.github.io/PESTO/reference/pesto_obs_schema.md)
  :

  Build a grounded `obs_schema` descriptor for a manifest

- [`as_manifest()`](https://max578.github.io/PESTO/reference/as_manifest.md)
  :

  Convert a PESTO ensemble result into a `pesto_ensemble_manifest`

- [`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md)
  : Write a manifest to YAML + sidecar data files

- [`read_manifest()`](https://max578.github.io/PESTO/reference/read_manifest.md)
  : Read a manifest from YAML + sidecar data files

- [`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
  : Verify the integrity of a manifest

## Visualisation

Publication-grade plots for convergence, ensembles, identifiability,
surrogate diagnostics.

- [`plot_phi()`](https://max578.github.io/PESTO/reference/plot_phi.md) :
  Plot Objective Function (Phi) Convergence
- [`plot_ensemble()`](https://max578.github.io/PESTO/reference/plot_ensemble.md)
  : Plot Ensemble Parameter Distributions
- [`plot_identifiability()`](https://max578.github.io/PESTO/reference/plot_identifiability.md)
  : Plot Parameter Identifiability
- [`plot_surrogate_diagnostics()`](https://max578.github.io/PESTO/reference/plot_surrogate_diagnostics.md)
  : Plot Surrogate Diagnostics
