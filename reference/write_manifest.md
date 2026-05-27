# Write a manifest to YAML + sidecar data files

Serialises a `pesto_ensemble_manifest` as a YAML file at `file` plus
three data sidecars (`<basename>_params.<ext>`,
`<basename>_outputs.<ext>`, `<basename>_assim.<ext>`) in the same
directory. `<ext>` depends on `format`:

## Usage

``` r
write_manifest(manifest, ...)
```

## Arguments

- manifest:

  A `pesto_ensemble_manifest`.

- ...:

  Method-specific arguments. For `pesto_ensemble_manifest`: `file`
  (character path to the YAML output; parent directory must exist;
  sidecars are written next to it) and `format` (one of `"rds"`,
  `"both"`, `"csv_unverified"`).

## Value

Invisible character vector of the written paths (YAML + sidecars, in
write order).

## Details

- `"rds"` (default) — RDS sidecars only. IEEE 754 doubles round-trip
  bit-exactly;
  [`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
  recomputes the SHA-256 hash and confirms integrity.

- `"both"` — RDS sidecars **plus** parallel inspection CSVs
  (`<basename>_params_inspection.csv`, etc.). The hash is still bound to
  the RDS form; the CSVs are decorative only and are recorded in the
  YAML's `inspection_csv:` block.

- `"csv_unverified"` — CSV sidecars only. The hash is still recorded
  (computed from the in-memory binary representation) but
  [`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
  cannot recompute it from disk: CSV write-formatter precision loss (~1
  ULP at IEEE 754 epsilon) would falsely fail the check. The YAML
  carries `integrity: not_verifiable` so downstream tools can branch
  accordingly. Renamed from `"csv"` in PESTO 0.3.2 (post critical
  review): the old name was indistinguishable at a glance from the
  verifiable modes, which the review judged a footgun. Old YAMLs with
  `format: csv` continue to read back correctly.
