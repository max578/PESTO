# Surrogate-accelerated iterative ensemble smoother

## Overview

The surrogate-accelerated IES reduces the number of expensive
forward-model evaluations during iterative parameter estimation. It
embeds an online Gaussian process (GP) surrogate in the IES update loop
and uses uncertainty-driven switching to decide which ensemble members
need a full model evaluation and which can use the cheaper surrogate
prediction. This vignette demonstrates the workflow, states where it
pays off, and shows the savings on a real run.

### Regime of applicability

The GP surrogate has a bounded operating envelope, worth stating before
the demo. As a soft floor, use `n_train >= 5 * n_params`: below that
ratio the GP posterior variance rarely drops far enough for the
switching rule to fire, and the surrogate defers to the full model on
most realisations. The second axis is smoothness: smooth, near-linear
responses are well approximated by an RBF-kernel GP and yield large
savings, whereas rough or near-discontinuous responses (sharp wetting
fronts, threshold activations, regime switches) defeat the kernel and
savings collapse. The third is ensemble size: larger ensembles improve
the GP fit roughly as \mathcal{O}(n^{-1}) in posterior variance, and
fewer than about twenty realisations rarely produce a usable surrogate.
Outside this envelope, run pure IES via
[`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md):
[`surrogate_ensemble_update()`](https://max578.github.io/PESTO/reference/surrogate_ensemble_update.md)
reports near-zero savings rather than degrading the posterior, but you
pay the GP training cost for no return. The
[`check_surrogate_regime()`](https://max578.github.io/PESTO/reference/check_surrogate_regime.md)
helper flags an unfavourable regime before you commit.

## The test problem

A nonlinear forward model inspired by groundwater flow, with a linear
sensitivity and a nonlinear interaction term:

y_j = \eta_j + \alpha \sin(\eta_j)\\ e^{-0.1\\\lvert \eta_j \rvert} +
\varepsilon_j, \qquad \eta_j = \sum_k G\_{jk}\theta_k .

Here \eta_j is the linear predictor and \alpha controls the
nonlinearity: \alpha = 0 is linear, \alpha = 1 is strongly nonlinear.

``` r

library(PESTO)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(ggplot2)

forward_model <- function(theta, G, alpha = 0.3) {
  linear <- as.numeric(G %*% theta)
  nonlinear <- alpha * sin(linear) * exp(-0.1 * abs(linear))
  linear + nonlinear
}

set.seed(42)
n_par <- 20
n_obs <- 50
n_real <- 50

theta_true <- rnorm(n_par)
G <- matrix(rnorm(n_obs * n_par, sd = 1 / sqrt(n_par)), n_obs, n_par)
y_obs <- forward_model(theta_true, G, alpha = 0.3) + rnorm(n_obs, sd = 0.1)
```

## Step 1: generate the prior ensemble

``` r

par_ens <- matrix(rnorm(n_real * n_par, sd = 1.5), n_real, n_par)

# Run the forward model for each realisation
obs_ens <- t(apply(par_ens, 1, function(p) {
  forward_model(p, G, alpha = 0.3) + rnorm(n_obs, sd = 0.1)
}))

cat("Parameter ensemble:", nrow(par_ens), "realisations x", ncol(par_ens), "parameters\n")
#> Parameter ensemble: 50 realisations x 20 parameters
cat("Observation ensemble:", nrow(obs_ens), "realisations x", ncol(obs_ens), "observations\n")
#> Observation ensemble: 50 realisations x 50 observations
```

## Step 2: train the GP surrogate

PESTO trains a GP with an RBF kernel directly from the ensemble:

``` r

gp <- train_gp_surrogate(par_ens, obs_ens)

cat("Training samples:", gp$n_train, "\n")
#> Training samples: 50
cat("Parameters:", gp$n_par, "\n")
#> Parameters: 20
cat("Outputs:", gp$n_obs, "\n")
#> Outputs: 50
cat("Length scale:", round(gp$length_scale, 4), "\n")
#> Length scale: 9.2136
cat("Signal variance:", round(gp$signal_var, 4), "\n")
#> Signal variance: 2.5654
cat("Log marginal likelihood:", round(gp$log_marginal_likelihood, 1), "\n")
#> Log marginal likelihood: -3043.3
```

## Step 3: evaluate surrogate predictions

``` r

predictions <- predict_gp_surrogate(gp, par_ens)

cat("Prediction dimensions:", dim(predictions$mean), "\n")
#> Prediction dimensions: 50 50
cat("Mean prediction uncertainty:", round(mean(predictions$uncertainty), 6), "\n")
#> Mean prediction uncertainty: 0.009998
cat("Max prediction uncertainty:", round(max(predictions$uncertainty), 6), "\n")
#> Max prediction uncertainty: 0.009999

# Check prediction accuracy
pred_error <- sqrt(mean((predictions$mean - obs_ens)^2))
cat("Prediction RMSE:", round(pred_error, 6), "\n")
#> Prediction RMSE: 0.000119
```

The GP achieves near-zero prediction error because it is interpolating
within its training set – this is by design, not an artefact.

## Step 4: surrogate-assisted IES update

[`surrogate_ensemble_update()`](https://max578.github.io/PESTO/reference/surrogate_ensemble_update.md)
performs the full algorithm:

1.  Train a GP from the current ensemble.
2.  Predict model outputs for all realisations.
3.  Classify by uncertainty: high -\> full model, low -\> surrogate.
4.  Apply the control-variate bias correction.
5.  Compute the IES update on the blended ensemble.

``` r

weights <- rep(1 / 0.1, n_obs)    # 1 / obs_noise
parcov_inv <- rep(1 / 1.5^2, n_par)  # 1 / prior_var

result <- surrogate_ensemble_update(
  par_ensemble = par_ens,
  obs_ensemble = obs_ens,
  obs_target = y_obs,
  weights = weights,
  parcov_inv = parcov_inv,
  cur_lam = 1.0,
  uncertainty_threshold = 0.1
)

cat("Full model runs needed:", result$n_model_runs, "\n")
#> Full model runs needed: 0
cat("Surrogate predictions:", result$n_surrogate_runs, "\n")
#> Surrogate predictions: 50
cat("Total realisations:", result$n_total, "\n")
#> Total realisations: 50
cat("Model savings:", sprintf("%.0f%%", result$savings_pct), "\n")
#> Model savings: 100%
cat("Mean GP uncertainty:", round(result$mean_uncertainty, 6), "\n")
#> Mean GP uncertainty: 0.009998
cat("GP length scale:", round(result$gp_length_scale, 4), "\n")
#> GP length scale: 9.2136
```

This twenty-parameter problem with a fifty-member ensemble sits *below*
the `n_train >= 5 * n_params` floor, so the GP variance rarely clears
the switching threshold and the savings here are modest by design. The
next section places the same machinery in a favourable regime, where it
pays off.

## Comparison: standard vs surrogate

``` r

# Standard IES update (for comparison)
par_mean <- colMeans(par_ens)
obs_mean <- colMeans(obs_ens)
par_diff <- t(par_ens) - par_mean
obs_diff <- t(obs_ens) - obs_mean
# Sign: ensemble_solution() expects sim - obs (see ?ensemble_solution).
obs_resid <- t(obs_ens) - matrix(rep(y_obs, n_real), n_obs, n_real)
par_resid <- par_diff
Am <- matrix(rnorm(n_par * (n_real - 1)), n_par, n_real - 1)

standard_upgrade <- ensemble_solution(
  par_diff, obs_diff, obs_resid, par_resid,
  weights, parcov_inv, Am, cur_lam = 1.0
)

surrogate_upgrade <- result$upgrade

# Compare upgrades
upgrade_diff <- sqrt(mean((standard_upgrade - surrogate_upgrade)^2))
upgrade_corr <- cor(as.numeric(standard_upgrade), as.numeric(surrogate_upgrade))

cat("Upgrade RMSE difference:", round(upgrade_diff, 6), "\n")
#> Upgrade RMSE difference: 0.000118
cat("Upgrade correlation:", round(upgrade_corr, 6), "\n")
#> Upgrade correlation: 1
```

## When the surrogate pays off

The saving is the fraction of realisations the surrogate handles, 1 -
n\_{\text{model}}/n\_{\text{total}}. Because the GP is trained on the
current ensemble and predicts at those same points, in a smooth regime
it is confident enough to handle essentially the whole ensemble at the
default threshold – so the saving is high *by construction*. What
decides whether the surrogate is worth using is therefore not how many
realisations it handles but whether its answer stays accurate, which the
*Effect of nonlinearity* section measures next.

``` r

set.seed(7)
np_f  <- 4L
no_f  <- 12L
nr_f  <- 120L
G_f   <- matrix(rnorm(no_f * np_f, sd = 1 / sqrt(np_f)), no_f, np_f)
fwd_f <- function(theta) as.numeric(G_f %*% theta)   # smooth, linear response
th_f  <- rnorm(np_f)
y_f   <- fwd_f(th_f) + rnorm(no_f, sd = 0.05)
pe_f  <- matrix(rnorm(nr_f * np_f, sd = 1.0), nr_f, np_f)
oe_f  <- t(apply(pe_f, 1L, function(p) fwd_f(p) + rnorm(no_f, sd = 0.05)))

r_f    <- surrogate_ensemble_update(pe_f, oe_f, y_f, rep(1 / 0.05, no_f),
                                    rep(1, np_f), uncertainty_threshold = 0.1)
post_f <- pe_f + r_f$upgrade
cat(sprintf(
  "Smooth 4-parameter problem: %.0f%% of evaluations handled by the surrogate; posterior RMSE %.3f\n",
  r_f$savings_pct, sqrt(mean((colMeans(post_f) - th_f)^2))))
#> Smooth 4-parameter problem: 100% of evaluations handled by the surrogate; posterior RMSE 0.014
```

The corollary matters as much: in the curse-of-dimensionality regime
(parameters outnumbering the ensemble) the surrogate still predicts the
whole ensemble, but its answer degrades – high savings with a poor
posterior. That is why
[`check_surrogate_regime()`](https://max578.github.io/PESTO/reference/check_surrogate_regime.md)
warns before you commit, and why the accuracy check below is the one
that counts.

## Effect of nonlinearity

How does surrogate accuracy vary with model nonlinearity?

``` r

results_dt <- rbindlist(lapply(seq(0, 1, by = 0.1), function(alpha) {
  set.seed(42)
  y_a <- forward_model(theta_true, G, alpha) + rnorm(n_obs, sd = 0.1)
  obs_a <- t(apply(par_ens, 1, function(p) {
    forward_model(p, G, alpha) + rnorm(n_obs, sd = 0.1)
  }))

  pm <- colMeans(par_ens); om <- colMeans(obs_a)
  pd <- t(par_ens) - pm; od <- t(obs_a) - om
  or_ <- t(obs_a) - matrix(rep(y_a, n_real), n_obs, n_real)
  pr <- pd
  Am_a <- matrix(rnorm(n_par * (n_real - 1)), n_par, n_real - 1)

  upg_std <- ensemble_solution(pd, od, or_, pr, weights, parcov_inv, Am_a, 1.0)
  par_std <- par_ens + upg_std
  rmse_std <- sqrt(mean((colMeans(par_std) - theta_true)^2))

  surr <- surrogate_ensemble_update(par_ens, obs_a, y_a, weights, parcov_inv,
                                     uncertainty_threshold = 0.1)
  par_surr <- par_ens + surr$upgrade
  rmse_surr <- sqrt(mean((colMeans(par_surr) - theta_true)^2))

  data.table(alpha = alpha, rmse_std = rmse_std, rmse_surr = rmse_surr,
             ratio = rmse_surr / rmse_std, savings = surr$savings_pct)
}))

ggplot(results_dt, aes(x = alpha)) +
  geom_line(aes(y = ratio), colour = "#009E73", linewidth = 1.2) +
  geom_point(aes(y = ratio), colour = "#009E73", size = 3) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0.95, 1.05)) +
  labs(title = "Surrogate accuracy across nonlinearity levels",
       x = expression(paste("Nonlinearity strength ", alpha)),
       y = "RMSE ratio (surrogate / standard)",
       caption = "Ratio = 1.0 means identical posterior quality") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))
```

![Surrogate RMSE ratio across nonlinearity
levels.](surrogate-ies_files/figure-html/nonlinearity-sweep-1.png)

Surrogate RMSE ratio across nonlinearity levels.

## Random Fourier features for large ensembles

The exact GP has O(n^3) training cost. For large ensembles (n \> 300),
the random Fourier feature (RFF) approximation scales as O(nD^2), where
D is the number of random features:

``` r

# Train RFF surrogate (D = 200 features)
rff <- train_rff_surrogate(par_ens, obs_ens, n_features = 200L)
cat("Training MSE:", round(rff$train_mse, 6), "\n")
#> Training MSE: 0

# Predict
rff_pred <- predict_rff_surrogate(rff, par_ens)
rff_error <- sqrt(mean((rff_pred$mean - obs_ens)^2))
cat("RFF prediction RMSE:", round(rff_error, 4), "\n")
#> RFF prediction RMSE: 5e-04

# Compare GP vs RFF accuracy
gp_error <- sqrt(mean((predictions$mean - obs_ens)^2))
cat("Exact GP prediction RMSE:", round(gp_error, 6), "\n")
#> Exact GP prediction RMSE: 0.000119
cat("RFF / GP error ratio:", round(rff_error / max(gp_error, 1e-10), 1), "x\n")
#> RFF / GP error ratio: 4.5 x
```

## Scaling comparison

``` r

scaling_dt <- rbindlist(lapply(c(30, 50, 100, 200, 500), function(nr) {
  X <- matrix(rnorm(nr * 20), nr, 20)
  Y <- matrix(rnorm(nr * 50), nr, 50)

  t_gp <- system.time(train_gp_surrogate(X, Y))[["elapsed"]] * 1000
  t_rff <- system.time(train_rff_surrogate(X, Y, 200L))[["elapsed"]] * 1000

  data.table(n = nr, `Exact GP` = t_gp, `RFF (D=200)` = t_rff)
}))

scaling_long <- melt(scaling_dt, id.vars = "n",
                     variable.name = "method", value.name = "time_ms")

ggplot(scaling_long, aes(x = n, y = pmax(time_ms, 0.01), colour = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_y_log10() +
  scale_colour_manual(values = c("Exact GP" = "#D55E00", "RFF (D=200)" = "#0072B2")) +
  labs(title = "GP surrogate training time",
       x = "Ensemble size", y = "Training time (ms, log scale)",
       colour = "Method") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
```

![GP vs RFF training-time
scaling.](surrogate-ies_files/figure-html/scaling-1.png)

GP vs RFF training-time scaling.

## When does the surrogate help?

The surrogate is most effective when:

- the ensemble occupies a **compact region** of parameter space (GP
  interpolation is accurate);
- the forward model varies **smoothly** over the parameter range;
- the ensemble is **large enough** to train the GP (practically n \geq
  20).

It degrades gracefully when:

- the parameter space is very high-dimensional (n_p \gg 100);
- the model has **sharp discontinuities** or bifurcations;
- the ensemble spans a very **wide region** (uninformative prior).

In all cases, the uncertainty threshold \tau makes the algorithm fall
back to full model evaluations when the surrogate is unreliable.

## References

- Rasmussen, C. E. & Williams, C. K. I. (2006). *Gaussian Processes for
  Machine Learning*. MIT Press.
- Rahimi, A. & Recht, B. (2007). Random features for large-scale kernel
  machines. *Advances in Neural Information Processing Systems* 20.
- Chen, Y. & Oliver, D. S. (2013). Levenberg-Marquardt forms of the
  iterative ensemble smoother. *Computational Geosciences*, 17(4),
  689–703.
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*.
  Springer.
