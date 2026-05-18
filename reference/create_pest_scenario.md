# Create a PEST Scenario Programmatically

Builds a `pesto_pst` object from R data structures, without requiring an
existing .pst file.

## Usage

``` r
create_pest_scenario(
  parameters,
  observations,
  model_command,
  template_files = NULL,
  instruction_files = NULL,
  pestpp_options = list()
)
```

## Arguments

- parameters:

  data.table. Parameter definitions with columns: `parnme`, `partrans`,
  `parchglim`, `parval1`, `parlbnd`, `parubnd`, `pargp`.

- observations:

  data.table. Observation definitions with columns: `obsnme`, `obsval`,
  `weight`, `obgnme`.

- model_command:

  Character. Model command line(s).

- template_files:

  data.table. Template/model input file pairs with columns
  `template_file`, `model_file`.

- instruction_files:

  data.table. Instruction/model output file pairs with columns
  `instruction_file`, `model_file`.

- pestpp_options:

  Named list. PEST++ options.

## Value

A `pesto_pst` object.

## Examples

``` r
pars <- data.table::data.table(
  parnme = c("k1", "k2"),
  partrans = "log",
  parchglim = "factor",
  parval1 = c(1.0, 0.5),
  parlbnd = c(0.01, 0.001),
  parubnd = c(100, 50),
  pargp = "hydraulic"
)
obs <- data.table::data.table(
  obsnme = c("h1", "h2"),
  obsval = c(1.0, 2.0),
  weight = c(1.0, 1.0),
  obgnme = "head"
)
pst <- create_pest_scenario(
  parameters    = pars,
  observations  = obs,
  model_command = "python model.py"
)
inherits(pst, "pesto_pst")
#> [1] TRUE
pst$control_data$npar
#> [1] 2
```
