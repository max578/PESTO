# Run PEST++ IES (Iterative Ensemble Smoother)

Executes the pestpp-ies algorithm on a PEST control file. This is the
primary method for ensemble-based parameter estimation and uncertainty
quantification.

## Usage

``` r
pesto_ies(
  pst_file,
  exe = NULL,
  num_reals = 50,
  noptmax = NULL,
  lambda_scale_fac = c(0.1, 0.5, 1),
  ies_par_en = NULL,
  extra_args = list(),
  working_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- pst_file:

  Character. Path to the .pst control file.

- exe:

  Character. Path to the pestpp-ies executable. `NULL` (the default)
  resolves it from `PESTPP_IES_EXE_PATH`, then `PESTPP_BIN_DIR`, then
  the `PATH`. PESTO does not bundle PEST++; it drives an installation
  you supply.

- num_reals:

  Integer. Number of ensemble realisations (overrides the value in the
  .pst file).

- noptmax:

  Integer or `NULL`. Maximum number of iterations, written to field 1 of
  line 7 of the control file's `* control data` section. `NULL` (the
  default) leaves the control file's own value alone.

- lambda_scale_fac:

  Numeric vector. Lambda scaling factors.

- ies_par_en:

  Character. Path to existing parameter ensemble file.

- extra_args:

  Named list. Additional PEST++ options, written to the control file as
  `++key(value)` lines. Keys must be PestppOptions keys – PEST++ rejects
  one it does not recognise. An option the control file already sets is
  replaced, not duplicated.

- working_dir:

  Character. Working directory for the run. Defaults to the directory
  containing the .pst file.

- verbose:

  Logical. Print stdout/stderr from pestpp-ies.

## Value

A list of class `pesto_ies_result` containing:

- phi:

  data.table of objective function values per iteration

- par_ensemble:

  Final parameter ensemble (data.table)

- obs_ensemble:

  Final observation ensemble (data.table)

- exit_code:

  Integer exit code from pestpp-ies

- runtime_seconds:

  Total wall-clock runtime

## Examples

``` r
# \donttest{
if (nzchar(Sys.which("pestpp-ies"))) {
  pars <- data.table::data.table(
    parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
    parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
    parubnd = c(100, 50), pargp = "hydraulic"
  )
  obs <- data.table::data.table(
    obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
    weight = c(1.0, 1.0), obgnme = "head"
  )
  pst <- create_pest_scenario(pars, obs, model_command = "echo run")
  tf <- tempfile(fileext = ".pst")
  on.exit(unlink(tf), add = TRUE)
  write_pst(pst, tf)
  res <- pesto_ies(tf, num_reals = 3, noptmax = 1, verbose = FALSE)
  res$exit_code
}
# }
```
