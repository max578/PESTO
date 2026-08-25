# Read a manifest from YAML + sidecar data files

Inverse of
[`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md).
Reads the YAML, loads the three sidecar data files (paths resolved
relative to the YAML file), and reconstructs the
`pesto_ensemble_manifest` S7 object. The file extensions in the YAML's
`artefacts:` block determine the read path (`.rds` via `readRDS`, `.csv`
via [`utils::read.csv`](https://rdrr.io/r/utils/read.table.html)).

## Usage

``` r
read_manifest(file)
```

## Arguments

- file:

  Character. Path to the YAML manifest file.

## Value

A `pesto_ensemble_manifest`.

## Examples

``` r
m <- as_manifest(
  pesto_ies_callback(
    function(t) t %*% t(matrix(1, 2L, 2L)),
    matrix(rnorm(20L), 10L, 2L, dimnames = list(NULL, c("a", "b"))),
    c(o1 = 0.1, o2 = -0.2), obs_sd = 0.1, noptmax = 2L, verbose = FALSE
  ),
  seed = 1L
)
f <- tempfile(fileext = ".yaml")
write_manifest(m, f)
m2 <- read_manifest(f)
class(m2)
#> [1] "PESTO::pesto_ensemble_manifest" "S7_object"                     
```
