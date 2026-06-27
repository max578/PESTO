# Developer-side generator for the cross-tool benchmark golden.
#
# Distils the frozen, real outputs of the standalone PEST / PEST++ / PESTO
# reproducibility harness into a small, self-describing summary that is
# SHIPPED in inst/extdata/pestpp_cache/ and rendered by the comparison
# vignette. The shipped artefact carries scientific content plus version /
# date provenance only -- no developer filesystem paths.
#
# Run from the package root with:
#   Rscript tools/pestpp_benchmark/build_multitool_summary.R
#
# Source results (developer-side, .Rbuildignored harness; override with the
# PESTO_BENCHMARK_RESULTS environment variable):
#   <harness>/results/paper_v2_20260627/summary_by_tool.csv
#
# Provenance of the source run (verified from all_runs.csv tool_version):
#   PEST (classic)        18.25   (pest_predunc = GLM + PREDUNC; pest_nsmc = null-space MC)
#   PEST++ (pestpp-ies)   5.2.16
#   PESTO                 0.8.0.9000
#   20 fixed seeds per (tool x problem) cell; run 2026-06-27.
#
# Re-run policy: regenerate only when the upstream benchmark run is
# refreshed. The golden is frozen real data -- never hand-edit the numbers.

suppressMessages(library(data.table))

src_dir <- Sys.getenv(
  "PESTO_BENCHMARK_RESULTS",
  unset = file.path("..", "benchmark_pest_pestpp_pesto",
                    "results", "paper_v2_20260627")
)
csv_path <- file.path(src_dir, "summary_by_tool.csv")
if (!file.exists(csv_path)) {
  stop("Source summary not found: ", csv_path,
       "\nSet PESTO_BENCHMARK_RESULTS to the harness results directory.",
       call. = FALSE)
}

raw <- data.table::fread(csv_path)

# Friendly, public method labels (the shipped artefact must not leak the
# internal tool ids' provenance beyond the documented version table).
method_label <- c(
  pesto_callback  = "PESTO (native IES)",
  pesto_surrogate = "PESTO (surrogate-accelerated)",
  pest_predunc    = "PEST (classic, GLM + PREDUNC)",
  pest_nsmc       = "PEST (classic, null-space Monte Carlo)",
  pestpp_ies      = "PEST++ (pestpp-ies)"
)
method_family <- c(
  pesto_callback  = "PESTO",
  pesto_surrogate = "PESTO",
  pest_predunc    = "PEST",
  pest_nsmc       = "PEST",
  pestpp_ies      = "PEST++"
)

summary <- data.frame(
  tier        = raw$tier_id,
  tool        = raw$tool,
  method      = unname(method_label[raw$tool]),
  family      = unname(method_family[raw$tool]),
  n_seeds     = as.integer(raw$n_seeds),
  rmse_med    = raw$rmse_med,
  rmse_iqr_lo = raw$rmse_iqr_lo,
  rmse_iqr_hi = raw$rmse_iqr_hi,
  ci90_cov    = raw$ci90_cov_mean,
  ks_dist_med = raw$ks_dist_med,
  wallclock_s = raw$wallclock_med_sec,
  fwd_evals   = raw$fwd_evals_med,
  conv_rate   = raw$convergence_rate,
  stringsAsFactors = FALSE
)

problems <- data.frame(
  tier        = c("tier1", "tier2"),
  id          = c("linear_p20_n50", "logistic_exp_p10_n30"),
  regime      = c("linear, well-posed", "non-linear ODE"),
  n_params    = c(20L, 10L),
  n_obs       = c(50L, 30L),
  description = c(
    paste("Linear inverse problem y = X theta + noise, condition number ~10;",
          "20 parameters, 50 observations, Gaussian prior."),
    paste("Mildly non-linear predator-prey ODE (logistic growth + decaying",
          "input); 10 positive parameters inferred in log space, 30",
          "observations.")
  ),
  stringsAsFactors = FALSE
)

multitool_benchmark <- list(
  generated_on  = "2026-06-27",
  harness       = paste("standalone PEST / PEST++ / PESTO reproducibility",
                        "harness; fixed-seed Monte Carlo, 20 seeds per cell"),
  tool_versions = c(
    "PESTO"               = "0.8.0.9000",
    "PEST (classic)"      = "18.25",
    "PEST++ (pestpp-ies)" = "5.2.16"
  ),
  metrics_note  = paste(
    "All figures are medians over 20 fixed seeds. rmse is to the known",
    "truth in the inversion (log) space; ci90_cov is mean 90% credible",
    "interval coverage; ks_dist_med is the median Kolmogorov-Smirnov",
    "distance of the posterior to a reference; wallclock_s and fwd_evals",
    "are per single inversion. conv_rate is the harness's internal",
    "stopping flag and does not track solution accuracy."
  ),
  problems      = problems,
  summary       = summary
)

out_dir <- file.path("inst", "extdata", "pestpp_cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "multitool_benchmark_summary.rds")
saveRDS(multitool_benchmark, out_path)
cat("Wrote", out_path, "\n")
cat("Methods x problems:", nrow(summary), "rows\n")
print(summary[, c("tier", "method", "rmse_med", "wallclock_s")])
