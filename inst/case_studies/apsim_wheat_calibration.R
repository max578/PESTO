# =============================================================================
# APSIM Wheat calibration with PESTO -- reproducible case-study driver
# =============================================================================
#
# A synthetic-truth recovery experiment (an OSSE / "twin experiment"): set
# known-true APSIM parameters, simulate per-season wheat yields, add
# measurement noise, then recover the parameters with PESTO's iterative
# ensemble smoother driving APSIM in-process via apsim_callback(). Because the
# truth is known by construction, correctness is verifiable without any field
# data -- the honest way to demonstrate a calibration method.
#
# This is the runnable companion to the "Calibrating APSIM with PESTO"
# case-study vignette; the vignette ships the frozen output this script writes.
#
# REQUIREMENTS
#   * R packages: PESTO, apsimx (>= 2.7.0), data.table
#   * A working APSIM Next Gen installation. The script does not hardcode any
#     path; point it at your install with environment variables (a
#     non-standard / x64-under-Rosetta build needs all three):
#         PESTO_APSIM_EXE       full path to the APSIM `Models` executable
#         PESTO_APSIM_EXAMPLES  directory containing `Wheat.apsimx`
#         PESTO_DOTNET_ROOT     .NET runtime root (if APSIM needs an explicit one)
#     If unset, the script relies on apsimx's own auto-detection. If APSIM
#     cannot be found it stops with a clear message rather than failing oddly.
#
# USAGE
#   Rscript inst/case_studies/apsim_wheat_calibration.R [output_dir]
#
# Output: <output_dir>/apsim_wheat_calibration_result.rds (default: tempdir()).
# -----------------------------------------------------------------------------

suppressMessages({
  library(PESTO)
  library(data.table)
})
if (!requireNamespace("apsimx", quietly = TRUE)) {
  message("apsimx is not installed; install it to run this case-study. Skipping.")
  quit(save = "no", status = 0L)
}

args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1L) args[[1]] else tempdir()
SEED     <- 20260628L
set.seed(SEED)

# ---- Configure + locate APSIM (no hardcoded paths) --------------------------
if (nzchar(Sys.getenv("PESTO_DOTNET_ROOT"))) {
  Sys.setenv(DOTNET_ROOT = Sys.getenv("PESTO_DOTNET_ROOT"))
}
apsim_exe <- Sys.getenv("PESTO_APSIM_EXE")
apsim_exd <- Sys.getenv("PESTO_APSIM_EXAMPLES")
if (nzchar(apsim_exe) && nzchar(apsim_exd)) {
  apsimx::apsimx_options(exe.path = apsim_exe, examples.path = apsim_exd,
                         warn.versions = FALSE, warn.find.apsimx = FALSE)
}
apsim_ok <- tryCatch(
  length(suppressWarnings(apsimx::apsim_version(which = "inuse",
                                                verbose = FALSE))) > 0L,
  error = function(e) FALSE
)
wheat_dir <- apsim_exd
wheat_src <- if (nzchar(wheat_dir)) file.path(wheat_dir, "Wheat.apsimx") else ""
if (!apsim_ok || !nzchar(wheat_src) || !file.exists(wheat_src)) {
  message("A working APSIM install and Wheat.apsimx were not found. ",
          "Set PESTO_APSIM_EXE / PESTO_APSIM_EXAMPLES. Skipping.")
  quit(save = "no", status = 0L)
}

# ---- Calibration target: 3 parameters with known truth ----------------------
# Paths resolved against the real Wheat.apsimx; apsim_callback's default writer
# edits each via apsimx::edit_apsimx(node = "Other", parm.path = <path>).
field <- ".Simulations.Simulation.Field"
spec <- list(
  CN2Bare     = list(path = paste0(field, ".Soil.SoilWater.CN2Bare"),
                     truth = 73,  lo = 50,  hi = 95,  unit = "curve number",
                     role  = "soil runoff (water limitation)"),
  FertN       = list(path = paste0(field, ".Fertilise at sowing.Amount"),
                     truth = 160, lo = 40,  hi = 260, unit = "kg N/ha",
                     role  = "fertiliser nitrogen (N limitation)"),
  SowingDepth = list(path = paste0(field, ".Sow using a variable rule.SowingDepth"),
                     truth = 30,  lo = 10,  hi = 90,  unit = "mm",
                     role  = "sowing depth (weak lever)")
)
pnm       <- names(spec)
param_map <- lapply(spec, `[[`, "path")
extract_yield <- function(sim) as.numeric(sim$Yield)

# ---- Synthetic truth: run APSIM at the true parameters -----------------------
truth_dir <- tempfile("truth_"); dir.create(truth_dir)
file.copy(wheat_src, truth_dir)
for (p in pnm) {
  apsimx::edit_apsimx("Wheat.apsimx", src.dir = truth_dir, wrt.dir = truth_dir,
                      node = "Other", parm.path = spec[[p]]$path,
                      value = spec[[p]]$truth, overwrite = TRUE, verbose = FALSE)
}
truth_sim <- apsimx::apsimx("Wheat.apsimx", src.dir = truth_dir, value = "report")
y_truth     <- extract_yield(truth_sim)
nseason     <- length(y_truth)
season_year <- as.integer(format(as.Date(truth_sim$Date), "%Y"))

obs_sd_val <- 250                       # kg/ha measurement uncertainty
y_obs_all  <- y_truth + rnorm(nseason, sd = obs_sd_val)

# Calibrate on 65% of seasons; hold out the rest for out-of-sample validation.
cal_idx <- sort(sample(seq_len(nseason), size = floor(nseason * 0.65)))
val_idx <- setdiff(seq_len(nseason), cal_idx)

# ---- Forward model: APSIM in-process, parallel over realisations -------------
fm_raw <- apsim_callback(template = wheat_src, param_map = param_map,
                         output_extractor = extract_yield)
ncores <- parallel::detectCores()
fm <- pesto_forward_model(fm_raw, parallel = "multicore",
                          n_cores = max(1L, ncores - 2L))

# ---- Prior ensemble + IES ----------------------------------------------------
N     <- 40L
prior <- sapply(pnm, function(p) runif(N, spec[[p]]$lo, spec[[p]]$hi))
colnames(prior) <- pnm

obs_all <- y_obs_all; names(obs_all) <- paste0("s", seq_len(nseason))
sd_vec  <- rep(obs_sd_val, nseason)
sd_vec[val_idx] <- 1e6                  # held-out seasons excluded from the fit

message(sprintf("Running IES: N=%d, noptmax=6 (~%d APSIM runs)...", N, N * 7L))
fit <- pesto_ies_callback(forward_model = fm, prior_ensemble = prior,
                          obs = obs_all, obs_sd = sd_vec, noptmax = 6L,
                          verbose = TRUE)

# ---- Summaries: recovery + held-out validation -------------------------------
post <- as.matrix(fit$par_ensemble[, pnm, with = FALSE])
oe   <- as.matrix(fit$obs_ensemble[,
          grep("^s[0-9]+$", names(fit$obs_ensemble)), with = FALSE])
# Failed realisations (APSIM crashes on some parameter draws) appear as NA
# rows; drop them from the predictive quantiles.
pred_lo  <- apply(oe[, val_idx, drop = FALSE], 2L, quantile, 0.05, na.rm = TRUE)
pred_hi  <- apply(oe[, val_idx, drop = FALSE], 2L, quantile, 0.95, na.rm = TRUE)
val_cov  <- mean(y_truth[val_idx] >= pred_lo & y_truth[val_idx] <= pred_hi)

# ---- Persist frozen result (clean version string, no machine paths) ----------
result <- list(
  spec = spec, pnm = pnm, y_truth = y_truth, y_obs_all = y_obs_all,
  season_year = season_year, cal_idx = cal_idx, val_idx = val_idx,
  obs_sd = obs_sd_val, N = N, noptmax = 6L, prior = prior, fit = fit,
  post = post, obs_ensemble = oe, val_coverage_truth = val_cov,
  apsim_version = "APSIM 2026.5.8046.0",
  apsimx_version = as.character(packageVersion("apsimx")),
  pesto_version  = as.character(packageVersion("PESTO")),
  seed = SEED, generated_on = format(Sys.Date())
)
out_path <- file.path(out_dir, "apsim_wheat_calibration_result.rds")
saveRDS(result, out_path)
message("Wrote ", out_path)
