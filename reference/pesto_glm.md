# Run PEST++ GLM (Gauss-Levenberg-Marquardt)

Executes the pestpp-glm algorithm for deterministic parameter
estimation.

## Usage

``` r
pesto_glm(
  pst_file,
  exe = NULL,
  noptmax = 20,
  extra_args = list(),
  working_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- pst_file:

  Character. Path to the .pst control file.

- exe:

  Character. Path to pestpp-glm executable.

- noptmax:

  Integer. Maximum number of iterations.

- extra_args:

  Named list. Additional PEST++ arguments.

- working_dir:

  Character. Working directory.

- verbose:

  Logical. Print output.

## Value

A list of class `pesto_glm_result`.

## Examples

``` r
# \donttest{
if (nzchar(Sys.which("pestpp-glm"))) {
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
  res <- pesto_glm(tf, noptmax = 1, verbose = FALSE)
  res$exit_code
}
# }
```
