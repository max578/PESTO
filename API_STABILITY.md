# API stability policy

`PESTO` follows a published two-phase policy. The phase boundary is the
1.0.0 release.

## Pre-1.0 (current)

Versions `0.x.y` follow **additive-by-intent** evolution. Every minor
release (`0.x → 0.(x+1)`) adds new exports without breaking existing
signatures; patch releases fix bugs and refresh documentation. Track
record so far:

- `0.1.0` — first release (core IES / GLM solvers, SVD backends,
  surrogates, PEST++ integration, visualisation).
- `0.1.1` — sole-authorship consolidation; no API change.
- `0.2.0` — added
  [`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
  (in-process IES driver) and
  [`apsim_callback()`](https://AAGI-AUS.github.io/PESTO/reference/apsim_callback.md)
  (apsimx adapter). **No breaking changes.**
- `0.3.0` — added `pesto_ensemble_manifest` (S7 cross-package contract)
  plus
  [`as_manifest()`](https://AAGI-AUS.github.io/PESTO/reference/as_manifest.md),
  [`write_manifest()`](https://AAGI-AUS.github.io/PESTO/reference/write_manifest.md),
  [`read_manifest()`](https://AAGI-AUS.github.io/PESTO/reference/read_manifest.md),
  [`verify_manifest()`](https://AAGI-AUS.github.io/PESTO/reference/verify_manifest.md).
  **No breaking changes.**
- `0.3.1` — added `format = c("rds","both","csv")` to
  [`write_manifest()`](https://AAGI-AUS.github.io/PESTO/reference/write_manifest.md);
  backwards-compatible YAML.
- `0.3.2` — renamed `format = "csv"` to `"csv_unverified"` at every call
  site; legacy `"csv"` still accepted with a deprecation warning, so
  0.3.1 manifests round-trip cleanly.
- `0.3.3` — Makevars FLIBS portability fix; no API change.

The policy is “additive-by-intent” rather than “frozen”: pre-1.0
reserves the right to break an existing signature when a design flaw
surfaces, but every such change must be:

1.  Listed under a `## Breaking changes` heading in `NEWS.md`, first.
2.  Justified in the release notes.
3.  Where feasible, accompanied by a temporary back-compat shim (the
    v0.3.2 `csv` → `csv_unverified` rename is the template: legacy
    spelling still works with a deprecation warning).

If you depend on `PESTO` pre-1.0, pin the version in `renv.lock` or
`DESCRIPTION` (`Imports: PESTO (>= 0.3.3)`).

## 1.0 and after

From `1.0.0`, `PESTO` adopts **strict frozen-API additive evolution** —
the same policy as `glmnet`, `mgcv`, and other long-lived solver-style
packages:

- Signatures of exported functions never change in a
  backwards-incompatible way across major versions.
- New capability arrives via new entry points (new exports), never via
  changes to existing ones.
- If an existing function genuinely needs to be retired, it is marked
  with
  [`lifecycle::deprecate_warn()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html),
  kept for ≥ 2 minor versions, then promoted to
  [`lifecycle::deprecate_stop()`](https://lifecycle.r-lib.org/reference/deprecate_soft.html)
  for ≥ 1 more minor version, then removed in the next major release.
  Successors are signposted in the deprecation message.

This policy is chosen because `PESTO` is invoked as the IES / surrogate
backend from pipeline code that the maintainer cannot edit — kernR
consumes the `pesto_ensemble_manifest` S7 contract; APSIM ensemble loops
via
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
are written once and run for the life of a calibration campaign. Silent
breakage across versions would defeat that contract.

## Versioning and tags

`PESTO` uses [Semantic Versioning 2.0.0](https://semver.org/). Every
release is git-tagged `vX.Y.Z` on `main`; the tag is annotated and
carries the `NEWS.md` entry as its message.

## Reporting an unintended break

If you discover that a `PESTO` release has silently broken your
pipeline, open an issue at <https://github.com/AAGI-AUS/PESTO/issues>
with the version you upgraded from and to plus a small reproducible
example. Unintended breaks at any stage are treated as bugs.
