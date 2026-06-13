# Ensemble Manifests -- the Cross-Package Contract

## Why a manifest?

A PESTO ensemble run produces a moderately large pile of state:
parameter ensemble, simulated outputs, target observations, IES weights,
the RNG seed, the lambda schedule, runtime context. To let `kernR`,
`proxymix`, and the paper-writing pipeline consume that output **without
reaching into PESTO-specific list internals**, v0.3.0 introduces a
versioned S7 contract object – `pesto_ensemble_manifest` – plus YAML+CSV
persistence with SHA-256 integrity checking.

This is Year-1 §A5 of the UQ ag-stack roadmap. It is deliberately
infrastructure rather than methodology: the value is that every
downstream consumer reads the same object, with the same provenance
guarantees, on every run.

## Constructing a manifest from an IES run

``` r

library(PESTO)

set.seed(20260516)
npar <- 3L; nobs <- 6L; nreal <- 60L

G <- matrix(rnorm(nobs * npar), nobs, npar)
theta_true <- c(1.0, -0.5, 2.0)
y <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = 0.05)
names(y) <- paste0("o", seq_len(nobs))

prior <- matrix(rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, paste0("p", seq_len(npar))))

fit <- pesto_ies_callback(
  forward_model  = function(theta) theta %*% t(G),
  prior_ensemble = prior,
  obs            = y,
  obs_sd         = 0.05,
  noptmax        = 4L,
  verbose        = FALSE
)
```

Wrap the result:

``` r

m <- as_manifest(fit, seed = 20260516L,
                 apsim_version = NA_character_)
print(m)
#> <pesto_ensemble_manifest> schema 1.1.0
#>   run_id        : ies_callback_20260613_231746_a0e19c00
#>   method        : ies_callback  (noptmax=4)
#>   ensemble      : 60 realisations x 3 parameters | 6 observations
#>   failure rate  : 0.00%
#>   pesto version : 0.7.0  apsim: NA
#>   timestamp     : 2026-06-13T23:17:46+0000
#>   data hash     : sha256:e7b630ad06429b528fa6a57a4973894eb9bf2709a6a3ffeb01366ede48b78ed6
```

Slots are reachable via the standard S7 `@` accessor:

``` r

m@run_id
#> [1] "ies_callback_20260613_231746_a0e19c00"
m@data_hash
#> [1] "sha256:e7b630ad06429b528fa6a57a4973894eb9bf2709a6a3ffeb01366ede48b78ed6"
m@noptmax
#> [1] 4
head(m@params)
#>   real_name       p1         p2       p3
#> 1    real_1 1.029000 -0.4654747 2.011874
#> 2    real_2 1.026500 -0.4658594 2.013726
#> 3    real_3 1.032378 -0.4622419 2.014441
#> 4    real_4 1.031622 -0.4617477 2.015089
#> 5    real_5 1.034491 -0.4633193 2.014088
#> 6    real_6 1.034192 -0.4611789 2.017081
```

## Writing, reading, and verifying

[`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md)
emits the YAML plus three RDS sidecars (`*_params.rds`, `*_outputs.rds`,
`*_assim.rds`). RDS is used in preference to CSV so IEEE 754 doubles
round-trip bit-exactly – the SHA-256 integrity check would otherwise
trip on formatter precision loss:

``` r

dir <- tempfile("pesto_manifest_")
dir.create(dir)
paths <- write_manifest(m, file.path(dir, "wagga_2026_run01.yaml"))
basename(paths)
#> [1] "wagga_2026_run01.yaml"        "wagga_2026_run01_params.rds" 
#> [3] "wagga_2026_run01_outputs.rds" "wagga_2026_run01_assim.rds"
```

Read it back and confirm the integrity hash matches:

``` r

m2 <- read_manifest(file.path(dir, "wagga_2026_run01.yaml"))
verify_manifest(m2)$ok
#> [1] TRUE
```

A peek at what the YAML actually looks like (truncated):

``` r

cat(paste(readLines(file.path(dir, "wagga_2026_run01.yaml"))[1:14],
          collapse = "\n"))
#> schema_version: 1.1.0
#> run_id: ies_callback_20260613_231746_a0e19c00
#> data_hash: sha256:e7b630ad06429b528fa6a57a4973894eb9bf2709a6a3ffeb01366ede48b78ed6
#> format: rds
#> integrity: verifiable
#> obs_schema: ~
#> artefacts:
#>   params: wagga_2026_run01_params.rds
#>   outputs: wagga_2026_run01_outputs.rds
#>   assim: wagga_2026_run01_assim.rds
#> seed: 20260516
#> fidelity: ~
#> apsim_version: ~
#> pesto_version: 0.7.0
```

### Inspection CSVs (optional)

For workflows that need a human-readable view of the ensemble – quick
scans in Excel, ad-hoc exports for a domain collaborator, copy-paste
into a Quarto narrative –
[`write_manifest()`](https://max578.github.io/PESTO/reference/write_manifest.md)
accepts a `format = "both"` argument that emits CSV sidecars
**alongside** the RDS. The hash stays bound to the RDS (bit-exact
integrity preserved), and the YAML records the inspection paths under a
separate `inspection_csv:` block so consumers can ignore them:

``` r

dir2 <- tempfile("pesto_manifest_csv_"); dir.create(dir2)
paths2 <- write_manifest(m, file.path(dir2, "wagga_2026_run01.yaml"),
                         format = "both")
basename(paths2)
#> [1] "wagga_2026_run01.yaml"                  
#> [2] "wagga_2026_run01_params.rds"            
#> [3] "wagga_2026_run01_outputs.rds"           
#> [4] "wagga_2026_run01_assim.rds"             
#> [5] "wagga_2026_run01_params_inspection.csv" 
#> [6] "wagga_2026_run01_outputs_inspection.csv"
#> [7] "wagga_2026_run01_assim_inspection.csv"
```

A `format = "csv_unverified"` mode (renamed from `"csv"` in PESTO 0.3.2)
is available for one-way exports to non-R analysts where round-trip
integrity is not required. The YAML carries `integrity: not_verifiable`
so downstream tools can branch on the weaker contract. The
[`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
function returns `ok = NA` rather than a spurious `FALSE` because CSV
formatter precision loss (~1 ULP at IEEE 754 epsilon) is enough to flip
the hash:

``` r

dir3 <- tempfile("pesto_manifest_unverified_"); dir.create(dir3)
write_manifest(m, file.path(dir3, "snapshot.yaml"),
               format = "csv_unverified")
m_csv <- read_manifest(file.path(dir3, "snapshot.yaml"))
verify_manifest(m_csv)$ok      # NA -- see $message for why
#> [1] NA
```

If you need both human inspection AND verifiable integrity, use
`format = "both"`. The `csv_unverified` mode is deliberately named to
flag the weaker contract at every call-site. It is for export, not for
storage you intend to re-load and trust.

## Tamper-detection

If a downstream tool silently re-saves the outputs sidecar (accidental
edit, partial overwrite) the hash will not match:

``` r

out_rds <- file.path(dir, "wagga_2026_run01_outputs.rds")
df <- readRDS(out_rds)
df[1, 2] <- df[1, 2] + 1e-3       # perturb one cell
saveRDS(df, out_rds, version = 3L)

m3 <- read_manifest(file.path(dir, "wagga_2026_run01.yaml"))
v  <- verify_manifest(m3)
v$ok
#> [1] FALSE
```

[`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
returns the stored vs recomputed hashes so a downstream consumer can
fail fast and report the divergence cleanly.

## Cross-package contract

The manifest is the single object that `kernR` (HSIC identifiability,
DR-DATE counterfactuals, MMD posterior-predictive checks) and `proxymix`
(GMM density-ratio bridges) consume. By dispatching on the S7 class,
those packages never see PESTO-internal list shapes. They read
`m@params`, `m@outputs`, `m@weights`, `m@obs_target` through the
contract and let
[`verify_manifest()`](https://max578.github.io/PESTO/reference/verify_manifest.md)
gate on integrity.

The companion jstyle `outputs_manifest.yaml` (the project-level artefact
index) will reference per-run manifests like this one by relative path +
the same `data_hash`, so a project’s manifest tree is end-to-end
hash-verifiable.

## Reproducibility

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] PESTO_0.7.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3        cli_3.6.6          knitr_1.51         rlang_1.2.0       
#>  [5] xfun_0.58          otel_0.2.0         S7_0.2.2           textshaping_1.0.5 
#>  [9] jsonlite_2.0.0     data.table_1.18.4  glue_1.8.1         htmltools_0.5.9   
#> [13] ragg_1.5.2         sass_0.4.10        scales_1.4.0       rmarkdown_2.31    
#> [17] grid_4.6.0         evaluate_1.0.5     jquerylib_0.1.4    fastmap_1.2.0     
#> [21] yaml_2.3.12        lifecycle_1.0.5    compiler_4.6.0     RColorBrewer_1.1-3
#> [25] fs_2.1.0           Rcpp_1.1.1-1.1     farver_2.1.2       systemfonts_1.3.2 
#> [29] digest_0.6.39      R6_2.6.1           bslib_0.11.0       gtable_0.3.6      
#> [33] tools_4.6.0        pkgdown_2.2.0      ggplot2_4.0.3      cachem_1.1.0      
#> [37] desc_1.4.3
```
