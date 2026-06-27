# PESTO test suite

The suite runs under [`testthat`](https://testthat.r-lib.org/) (third edition)
and is organised so that the regular run is fast and self-contained, while a
heavier *extended* battery is available behind a single flag. It is structured
to address the rOpenSci statistical-software **General Standards** for testing
(`srr` codes `G5.*`); the standard-by-standard compliance notes live next to the
tests that satisfy them (see the `@srrstats` blocks in `tests/testthat/`).

## Running the tests

```r
# regular suite (fast; what CRAN and the default CI job run)
devtools::test()

# full statistical battery (adds many-seed recovery + finer scaling sweeps)
Sys.setenv(PESTO_EXTENDED_TESTS = "true")
devtools::test()
```

## Extended tests (`PESTO_EXTENDED_TESTS`)

Long-running members of the battery are gated by `skip_if_not_extended()`
(defined in `tests/testthat/helper-srr.R`), which reads the
`PESTO_EXTENDED_TESTS` environment variable. They run under the *same*
`testthat` framework as the regular tests and are simply skipped when the flag
is unset. To run them in GitHub Actions, add to the workflow `env` block:

```yaml
env:
  PESTO_EXTENDED_TESTS: true
```

The extended tests are:

| Test | What it adds | Approx. runtime |
|------|--------------|-----------------|
| `test-recovery-scaling.R` — "many seeds" | parameter recovery over 8 seeds (vs 3 in the regular run) | < 5 s |

No test downloads external data: every dataset is simulated inline from a fixed
`set.seed()`, so the extended tests need no network access and no large assets.

## Optional-dependency and platform skips

Some tests exercise integrations that depend on optional packages or hardware.
They skip cleanly (with a diagnostic) when the requirement is absent:

| Requirement | Tests affected | Skip mechanism |
|-------------|----------------|----------------|
| `apsimx` package + APSIM binary | `test-apsim-callback.R`, the APSIM arm of `test-adapter-contract.R` | `skip_if_not_installed("apsimx")` |
| `deSolve` package | the `desolve` solver arm of `test-ode-templates.R` | `skip_if_not_installed("deSolve")` |
| `yaml` package | manifest round-trip in `test-ies-filter.R` | `skip_if_not_installed("yaml")` |
| >= 2 CPU cores, non-Windows | multicore arm of `test-forward-model.R` | `skip_on_os("windows")`, core check, `skip_on_cran()` |

## Fixtures

`tests/testthat/fixtures/pestpp_ies_tier1_golden.rds` is a self-contained golden
record used by the correctness test `test-correctness-analytic.R`: the posterior
mean produced by fixed-version **pestpp-ies 5.2.16** on the frozen `linear_p20_n50`
benchmark problem, together with the exact problem inputs (forward operator,
observations, prior ensemble, seeds) so PESTO can be re-run on the identical
problem and compared. Its `$provenance` element records the originating
benchmark run, tool versions, and seeds.
