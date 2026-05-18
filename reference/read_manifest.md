# Read a manifest from YAML + sidecar data files

Inverse of
[`write_manifest()`](https://AAGI-AUS.github.io/PESTO/reference/write_manifest.md).
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
