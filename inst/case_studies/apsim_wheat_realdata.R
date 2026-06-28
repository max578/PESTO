# =============================================================================
# APSIM Wheat calibration to REAL observed data with PESTO (case-study Part 2)
# =============================================================================
#
# Calibrates two physiological parameters of APSIM Wheat -- radiation-use
# efficiency (RUE) and the cultivar BasePhyllochron (phenology) -- to the real
# observed biomass in apsimx::obsWheat, using PESTO's iterative ensemble
# smoother. Unlike Part 1 (synthetic-truth recovery) there is no known truth;
# the checks are (a) goodness of fit, (b) out-of-sample prediction on held-out
# dates, and (c) agreement of PESTO's posterior with the independent point
# optimum from apsimx::optim_apsimx() (RUE ~ 1.5, BasePhyllochron ~ 87.6).
#
# The model (Wheat-opt-ex.apsimx), weather (Ames.met) and data (obsWheat) all
# ship with the apsimx package, so only a working APSIM binary is needed.
#
# REQUIREMENTS / USAGE: as for apsim_wheat_calibration.R. Point at an APSIM
# install with PESTO_APSIM_EXE (and PESTO_DOTNET_ROOT if needed); apsimx finds
# the bundled model/data itself.
#   Rscript inst/case_studies/apsim_wheat_realdata.R [output_dir]
# -----------------------------------------------------------------------------

suppressMessages({ library(PESTO); library(data.table) })
if (!requireNamespace("apsimx", quietly = TRUE)) {
  message("apsimx is not installed; skipping."); quit(save = "no", status = 0L)
}
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1]] else tempdir()
SEED    <- 20260628L; set.seed(SEED)

if (nzchar(Sys.getenv("PESTO_DOTNET_ROOT"))) {
  Sys.setenv(DOTNET_ROOT = Sys.getenv("PESTO_DOTNET_ROOT"))
}
if (nzchar(Sys.getenv("PESTO_APSIM_EXE"))) {
  apsimx::apsimx_options(exe.path = Sys.getenv("PESTO_APSIM_EXE"),
                         warn.versions = FALSE, warn.find.apsimx = FALSE)
}
apsim_ok <- tryCatch(
  length(suppressWarnings(apsimx::apsim_version(which = "inuse",
                                                verbose = FALSE))) > 0L,
  error = function(e) FALSE)
if (!apsim_ok) {
  message("No working APSIM install found (set PESTO_APSIM_EXE). Skipping.")
  quit(save = "no", status = 0L)
}

# ---- Real observed data + bundled model (all from apsimx) --------------------
utils::data("obsWheat", package = "apsimx")
obs_dates <- as.Date(obsWheat$Date)
y_obs     <- obsWheat$Wheat.AboveGround.Wt          # measured biomass (g/m2)
nobs      <- length(y_obs)

ed <- system.file("extdata", package = "apsimx")
wf <- "Wheat-opt-ex.apsimx"
wd <- tempfile("apsim_rd_"); dir.create(wd)
file.copy(file.path(ed, "Ames.met"), wd)            # shared weather

spec <- list(
  RUE         = list(path = "Wheat.Leaf.Photosynthesis.RUE.FixedValue",
                     lo = 0.8, hi = 2.2, default = 1.2, apsimx_opt = 1.50,
                     unit = "g/MJ", role = "radiation use efficiency"),
  Phyllochron = list(path = "Wheat.Cultivars.USA.Yecora.BasePhyllochron",
                     lo = 60, hi = 160, default = 120, apsimx_opt = 87.6,
                     unit = "degC.day", role = "phenology (phyllochron)")
)
pnm <- names(spec); param_map <- lapply(spec, `[[`, "path")

# RUE and BasePhyllochron live in a Replacements node, so they need
# edit_apsimx_replacement() (node.string = parent path, parm = leaf), not the
# default node="Other" writer. Supply it via apsim_callback's param_writer hook.
writer_repl <- function(file, src.dir, node, value) {
  toks <- strsplit(node, ".", fixed = TRUE)[[1]]
  apsimx::edit_apsimx_replacement(
    file, src.dir = src.dir, wrt.dir = src.dir,
    node.string = paste(utils::head(toks, -1), collapse = "."),
    parm = utils::tail(toks, 1), value = value,
    overwrite = TRUE, verbose = FALSE)
}
# Biomass at the observation dates (the model reports daily).
extract_bio <- function(sim) {
  s <- as.data.table(sim); s[, D := as.Date(Date)]
  bcol <- grep("AboveGround.Wt", names(s), value = TRUE)[1]
  s[[bcol]][match(obs_dates, s$D)]
}

fm_raw <- apsim_callback(template = file.path(ed, wf), param_map = param_map,
                         output_extractor = extract_bio,
                         param_writer = writer_repl, workdir = wd)
ncores <- parallel::detectCores()
fm <- pesto_forward_model(fm_raw, parallel = "multicore",
                          n_cores = max(1L, ncores - 2L))

val_idx <- c(4L, 7L, 9L); cal_idx <- setdiff(seq_len(nobs), val_idx)
obs_all <- y_obs; names(obs_all) <- paste0("d", seq_len(nobs))
sd_vec  <- pmax(30, 0.10 * y_obs); sd_vec[val_idx] <- 1e6   # hold out from fit

N     <- 40L
prior <- sapply(pnm, function(p) runif(N, spec[[p]]$lo, spec[[p]]$hi))
colnames(prior) <- pnm

message("Running IES on real obsWheat data (N=40, noptmax=6)...")
fit <- pesto_ies_callback(fm, prior_ensemble = prior, obs = obs_all,
                          obs_sd = sd_vec, noptmax = 6L, verbose = TRUE)

post <- as.matrix(fit$par_ensemble[, pnm, with = FALSE])
oe   <- as.matrix(fit$obs_ensemble[,
          grep("^d[0-9]+$", names(fit$obs_ensemble)), with = FALSE])

result <- list(
  spec = spec, pnm = pnm, obs_dates = obs_dates, y_obs = y_obs,
  cal_idx = cal_idx, val_idx = val_idx, prior = prior, fit = fit, post = post,
  obs_ensemble = oe, obsWheat = obsWheat, sd_floor = 30, sd_frac = 0.10,
  apsim_version = "APSIM 2026.5.8046.0",
  apsimx_version = as.character(packageVersion("apsimx")),
  pesto_version  = as.character(packageVersion("PESTO")),
  seed = SEED, generated_on = format(Sys.Date())
)
out_path <- file.path(out_dir, "apsim_realdata_result.rds")
saveRDS(result, out_path)
message("Wrote ", out_path)
