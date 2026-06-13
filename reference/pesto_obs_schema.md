# Build a grounded `obs_schema` descriptor for a manifest

Constructs the optional `obs_schema` slot of a
[pesto_ensemble_manifest](https://max578.github.io/PESTO/reference/pesto_ensemble_manifest.md):
a machine-checkable statement of what each output and parameter column
*means* (its physical quantity and unit), so a downstream consumer can
verify two manifests are commensurable by name rather than trusting a
positional convention. This is the single Independent-Oracle anchor the
manifest contract previously lacked.

## Usage

``` r
pesto_obs_schema(outputs = NULL, params = NULL)
```

## Arguments

- outputs:

  Optional `data.frame` describing output columns, with columns `name`,
  `quantity`, `unit` (character) and optionally `verified_on` (`Date`),
  `oracle_kind`, `evidence_path`.

- params:

  Optional `data.frame` describing parameter columns, with columns
  `name`, `apsim_node`, `unit` and the same optional provenance columns.

## Value

A validated `list(outputs = , params = )` suitable for the `obs_schema`
argument of
[pesto_ensemble_manifest](https://max578.github.io/PESTO/reference/pesto_ensemble_manifest.md).
Either side may be omitted (`NULL`).

## Details

Each row optionally carries provenance per the orchestra provenance
vocabulary: `verified_on` (a `Date`, or `NA` for an unverified fact),
`oracle_kind`, and `evidence_path`. A column is *grounded* only when
`verified_on` is a non-`NA` date.

## See also

The manifest contract `inst/manifest_contract.md` and the orchestra
provenance vocabulary.

## Examples

``` r
# Describe two output columns and one parameter, grounding the yield
# column against a dated oracle and leaving the rest unverified.
schema <- pesto_obs_schema(
  outputs = data.frame(
    name        = c("yield", "biomass"),
    quantity    = c("grain yield", "above-ground biomass"),
    unit        = c("kg/ha", "kg/ha"),
    verified_on = c(as.Date("2026-06-01"), as.Date(NA)),
    stringsAsFactors = FALSE
  ),
  params = data.frame(
    name      = "rue",
    apsim_node = "Wheat.Leaf.Photosynthesis.RUE.FixedValue",
    unit      = "g/MJ",
    stringsAsFactors = FALSE
  )
)
str(schema, max.level = 2L)
#> List of 2
#>  $ outputs:'data.frame': 2 obs. of  6 variables:
#>   ..$ name         : chr [1:2] "yield" "biomass"
#>   ..$ quantity     : chr [1:2] "grain yield" "above-ground biomass"
#>   ..$ unit         : chr [1:2] "kg/ha" "kg/ha"
#>   ..$ verified_on  : Date[1:2], format: "2026-06-01" NA
#>   ..$ oracle_kind  : chr [1:2] NA NA
#>   ..$ evidence_path: chr [1:2] NA NA
#>  $ params :'data.frame': 1 obs. of  6 variables:
#>   ..$ name         : chr "rue"
#>   ..$ apsim_node   : chr "Wheat.Leaf.Photosynthesis.RUE.FixedValue"
#>   ..$ unit         : chr "g/MJ"
#>   ..$ verified_on  : Date[1:1], format: NA
#>   ..$ oracle_kind  : chr NA
#>   ..$ evidence_path: chr NA
```
