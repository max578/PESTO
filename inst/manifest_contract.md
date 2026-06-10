# PESTO ensemble-manifest contract (v1.1.0)

<!-- constellation-contract: machine-readable header; charter at flexyBayes_dev/CONSTELLATION.md -->
```yaml
contract: pesto_ensemble_manifest
id: C2
version: 1.1.0
kind: data
owner: PESTO
consumers: [kernR, flexyBayes]
conformance_test: "verify_manifest() integrity check + write/read round-trip + obs_schema correspondence (where both manifests carry it)"
charter: flexyBayes_dev/CONSTELLATION.md
```

The versioned, hashed, provenance-tracked object PESTO emits so independent
consumers can validate and generalise a calibrated ensemble **without depending
on PESTO's internals**. This is the data contract for the calibration -> validation
(kernR) and calibration -> generalisation (flexyBayes) edges of the constellation.

## What it is

`pesto_ensemble_manifest` is an S7 object carrying the calibrated parameter
ensemble (`params`), the simulated output ensemble (`outputs`), observation
weights (`weights`), the calibration target (`obs_target`), a SHA-256
`data_hash` over the payload, and provenance (`method`, `failure_rate`,
`pesto_version`, `schema_version`, `run_id`, timestamp). The authoritative schema
and accessors live in `R/manifest.R` -- this document is the *contract surface*,
not a second copy of the schema.

## Contract obligations (producer side)

- Treat the manifest properties as a **public, versioned contract**. Bump
  `schema_version` (currently `1.1.0`) on any breaking change; consumers branch
  on it. `1.1.0` added the optional `obs_schema` descriptor over `1.0.0`
  additively (a `1.0.0` manifest reads back as `1.1.0` with `obs_schema = NULL`).
- `format = "rds"` (default) is integrity-verifiable: `verify_manifest()`
  recomputes the hash from disk. `format = "csv_unverified"` records the hash but
  cannot recompute it (CSV precision loss) and carries `integrity: not_verifiable`
  so consumers refuse rather than silently trust.
- **Grounded semantics (`obs_schema`, recommended).** State each output /
  parameter column's physical quantity + unit via `pesto_obs_schema()` so column
  meaning is machine-checkable, not positional. Per the Independent Oracle
  Principle, this is the field that lets a consumer detect a wrong-but-agreed
  convention. Provenance per the orchestra provenance vocabulary
  (`ORCHESTRA_dev/integration/provenance_vocabulary.md`): a column is *grounded*
  only when its `verified_on` is a non-`NA` date.

## Conformance (consumer side)

A consumer honours C2 by: accepting either an in-memory manifest or a path
(`read_manifest()`), calling `verify_manifest()` before use, refusing on
`integrity: not_verifiable` or a high `failure_rate`, and threading
`data_hash` / `pesto_version` / `run_id` into its own provenance so the lineage
closes (e.g. `kernR::mmd_ppc()` traces back to the producing run). When **both**
the consumed manifests carry an `obs_schema`, the consumer must additionally
refuse to compare them if a shared column disagrees on unit or quantity
(correspondence, not just structure) — see `kernR::dr_date_scenario()`.

## Edges (see the charter's status matrix)

- **PESTO -> kernR (live).** `mmd_ppc()`, `dr_date_scenario()` consume the manifest.
- **PESTO -> flexyBayes (planned).** `fb_pesto()` via the `ensemble_source` registry.
