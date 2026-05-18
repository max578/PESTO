# apsimx Forward-Model Adapter for PESTO IES

Builds an in-process forward-model closure for
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
that drives APSIM (Next Gen `.apsimx` or Classic `.apsim`) through the
`apsimx` package without going via the `.pst`-file path. Used as Year-1
§D4 of the UQ ag-stack roadmap (`uq_ag_stack_roadmap_v0.md`).

## Usage

``` r
apsim_callback(
  template,
  param_map,
  output_extractor,
  workdir = tempfile("apsim_cb_"),
  param_writer = NULL,
  simulation_runner = NULL,
  verbose = FALSE
)
```

## Arguments

- template:

  Character. Path to a working `.apsimx` (Next Gen) or `.apsim`
  (Classic) template file. Per-realisation copies are made into
  `workdir`; the template itself is never modified.

- param_map:

  Named list. Names are the parameter columns expected in `theta`;
  values are character node paths understood by the appropriate
  `apsimx::edit_*` function (e.g. `"Manager.SowingRule.Rule.Population"`
  for Next Gen).

- output_extractor:

  Function. Takes the object returned by
  [`apsimx::apsimx()`](https://rdrr.io/pkg/apsimx/man/apsimx.html) /
  `apsim()` (typically a data.frame of report variables) and returns a
  length-`nobs` numeric vector. The first successful realisation defines
  `nobs`.

- workdir:

  Character. Per-run working directory. Created if it does not exist;
  not cleaned automatically (so failures are inspectable). Default: a
  fresh `tempfile("apsim_cb_")`.

- param_writer:

  Optional function with signature
  `function(file, src.dir, node, value)`. Overrides the default apsimx
  editor dispatch. Useful for unusual node paths or non-standard apsimx
  versions.

- simulation_runner:

  Optional function with signature `function(file, src.dir)` returning
  the simulation object. Overrides the default
  [`apsimx::apsimx()`](https://rdrr.io/pkg/apsimx/man/apsimx.html) /
  `apsim()` dispatch.

- verbose:

  Logical. Forward verbose flag into apsimx calls and print
  per-realisation status (default `FALSE`).

## Value

A closure of signature `function(theta) -> obs` suitable for the
`forward_model =` argument of
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md).

## Details

The returned closure has signature `function(theta) -> obs`, where
`theta` is an `nreal x npar` matrix with column names matching
`names(param_map)`, and `obs` is an `nreal x nobs` matrix. Each row is
produced by:

1.  Copying `template` into a fresh per-realisation file under
    `workdir`.

2.  For each `(par_name, node_path)` in `param_map`, calling the
    appropriate `apsimx::edit_*` function to write `theta[i, par_name]`
    to that node.

3.  Calling
    [`apsimx::apsimx()`](https://rdrr.io/pkg/apsimx/man/apsimx.html) (or
    `apsim()` for Classic) and passing the returned simulation object to
    `output_extractor()`, which must return a length-`nobs` numeric
    vector.

Per-realisation failures (APSIM crash, missing output, extractor error)
populate an `NA` row;
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
then carries that realisation forward unchanged or aborts, depending on
its `on_failure` setting.

## APSIM version compatibility

Tested against the `apsimx` package API as of CRAN 2.7.x. The exact
editor function differs between APSIM Next Gen (`edit_apsimx*`) and
APSIM Classic (`edit_apsim`); selection is by file extension of
`template`. If your installed `apsimx` version exposes a different
editor signature, supply `param_writer` to override the default
per-parameter writer.

## Concurrency

Phase-1 D4 runs realisations **serially**. Parallel execution via
`future` or `mirai` is a planned follow-up; `apsimx`'s thread-safety
under ensemble load has not been verified.

## See also

[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
for the IES driver;
[`pesto_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies.md)
for the classic `.pst`-file path.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires apsimx and a working APSIM installation
fm <- apsim_callback(
  template  = "wheat_wagga.apsimx",
  param_map = list(
    RUE       = "Wheat.Leaf.Photosynthesis.RUE.FixedValue",
    CN2       = "Soil.SoilWater.CN2Bare"
  ),
  output_extractor = function(sim) {
    # sim is a data.frame; extract end-of-season yield trajectory
    as.numeric(sim$Wheat.Grain.Total.Wt)
  }
)
prior <- matrix(c(runif(40, 1.0, 2.0), runif(40, 60, 90)),
                ncol = 2, dimnames = list(NULL, c("RUE", "CN2")))
fit <- pesto_ies_callback(
  forward_model  = fm,
  prior_ensemble = prior,
  obs            = c(y1 = 4500, y2 = 5200),
  obs_sd         = 200,
  noptmax        = 4
)
} # }
```
