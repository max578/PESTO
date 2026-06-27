# PESTO ensemble-manifest format specification (v1.1.0)

The `pesto_ensemble_manifest` is the versioned, hashed, self-describing record
PESTO emits so that **any** downstream tool can validate and reuse a calibrated
ensemble *without depending on PESTO's internals*. This document is the format
specification (the data contract); the authoritative schema and accessors live
in `R/manifest.R`.

## What it is

`pesto_ensemble_manifest` is an S7 object carrying the calibrated parameter
ensemble (`params`), the simulated output ensemble (`outputs`), observation
weights (`weights`), the calibration target (`obs_target`), a SHA-256
`data_hash` over the payload, and provenance (`method`, `failure_rate`,
`pesto_version`, `schema_version`, `run_id`, timestamp, and the optional
`fidelity` and `obs_schema`). It persists as a YAML index plus data sidecars
(`<basename>_params.*`, `<basename>_outputs.*`, `<basename>_assim.*`).

## Producer obligations (PESTO)

- Treat the manifest properties as a **public, versioned contract**. Bump
  `schema_version` (currently `1.1.0`) on any breaking change; consumers branch
  on it. `1.1.0` added the optional `obs_schema` descriptor over `1.0.0`
  additively (a `1.0.0` manifest reads back as `1.1.0` with `obs_schema = NULL`).
- `format = "rds"` (default) is integrity-verifiable: `verify_manifest()`
  recomputes the hash from disk. `format = "csv_unverified"` records the hash but
  cannot recompute it (CSV precision loss) and carries `integrity: not_verifiable`
  so consumers refuse rather than silently trust.
- **Grounded semantics (`obs_schema`, recommended).** State each output /
  parameter column's physical quantity and unit via `pesto_obs_schema()` so
  column meaning is machine-checkable rather than positional. Per the
  Independent Oracle Principle, this is the field that lets a consumer detect a
  wrong-but-agreed convention. A column is *grounded* only when its
  `verified_on` is a non-`NA` date.

## Consumer conformance (any downstream tool)

A conformant consumer:

1. accepts either an in-memory manifest or a path (`read_manifest()`);
2. calls `verify_manifest()` before use, and refuses on
   `integrity: not_verifiable` or an unacceptably high `failure_rate`;
3. reads the typed slots (`params`, `outputs`, `weights`, `obs_target`) and
   branches on `schema_version` for forward compatibility;
4. threads `data_hash` / `pesto_version` / `run_id` into its own provenance so
   the lineage closes back to the producing run;
5. when comparing two manifests that both carry an `obs_schema`, refuses the
   comparison if a shared column disagrees on unit or quantity (correspondence,
   not merely structure).

**Conformance test:** `verify_manifest()` integrity check, plus a write/read
round-trip, plus (where both manifests carry it) `obs_schema` correspondence.
