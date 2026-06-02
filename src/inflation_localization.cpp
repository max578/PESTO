// PESTO: Scalable Tools for Estimation of Parameters
// Finite-ensemble pathology countermeasures: spread diagnostics,
// covariance localisation (correlation-based + Gaspari-Cohn), and the
// explicit-gain localised ensemble update.
//
// These kernels address ensemble under-dispersion and spurious
// finite-sample correlations in the iterative ensemble smoother (IES).
// The standard GLM kernel ensemble_solution() works in the reduced
// observation-anomaly subspace and never materialises the explicit
// npar x nobs Kalman gain; state-space localisation needs that gain, so
// it lives here in ensemble_solution_localized() rather than as a flag on
// the SVD path.
//
// Licensed under GPL-3.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::DiagonalMatrix;

//' Spectral Spread Effective Sample Size of a Parameter Ensemble
//'
//' Diagnoses ensemble collapse (under-dispersion) by the participation
//' ratio of the parameter-anomaly covariance eigenspectrum. Given the
//' parameter anomalies \eqn{\Delta\Theta} (deviations from the ensemble
//' mean, \code{npar x nreal}), the anomaly covariance is
//' \eqn{C = \Delta\Theta \Delta\Theta^{T} / (N - 1)} with eigenvalues
//' \eqn{\lambda_i = s_i^2 / (N - 1)}, where \eqn{s_i} are the singular
//' values of \eqn{\Delta\Theta}.
//'
//' The spectral spread-ESS is the participation ratio
//' \deqn{\mathrm{ESS} = \frac{(\sum_i \lambda_i)^2}{\sum_i \lambda_i^2}
//'                    = \frac{(\sum_i s_i^2)^2}{\sum_i s_i^4},}
//' the effective number of directions carrying variance. It is bounded in
//' \eqn{[1, r_{\max}]} with \eqn{r_{\max} = \min(\mathrm{npar}, N - 1)}:
//' equal to \eqn{r_{\max}} when variance is spread isotropically across all
//' modes, and approaching 1 when the ensemble collapses onto a single
//' direction. Because the ratio is invariant to a global rescaling of the
//' anomalies, it isolates the *shape* of the collapse (directional
//' degeneracy) from its *magnitude*; magnitude is tracked separately by the
//' R-side spread-retention ratio.
//'
//' @param par_diff Matrix (npar x nreal). Parameter anomalies (deviations
//'   from the ensemble mean). At least 2 columns are required.
//' @return A list with components \code{ess} (the spectral spread-ESS),
//'   \code{r_max} (the maximum attainable value
//'   \eqn{\min(\mathrm{npar}, N - 1)}), and \code{ess_ratio}
//'   (\code{ess / r_max}, in \eqn{(0, 1]}).
//' @references
//' Bretherton, C.S., Widmann, M., Dymnikov, V.P., Wallace, J.M. & Blade, I.
//' (1999). The effective number of spatial degrees of freedom of a
//' time-varying field. \emph{Journal of Climate}, 12(7), 1990--2009.
//' @examples
//' set.seed(1L)
//' # Healthy isotropic spread -> ESS near r_max
//' good <- matrix(rnorm(6L * 40L), 6L, 40L)
//' ensemble_spread_ess(good)$ess_ratio
//' # Collapsed onto one direction -> ESS near 1
//' v <- rnorm(6L)
//' bad <- outer(v, rnorm(40L)) + matrix(rnorm(6L * 40L, sd = 1e-3), 6L, 40L)
//' ensemble_spread_ess(bad)$ess_ratio
//' @export
// [[Rcpp::export]]
Rcpp::List ensemble_spread_ess(const Eigen::MatrixXd& par_diff)
{
    const int npar  = par_diff.rows();
    const int nreal = par_diff.cols();
    if (nreal < 2) {
        Rcpp::stop("`par_diff` must have at least 2 columns (realisations).");
    }

    // Singular values only (no U / V): cheaper than the thin factors.
    Eigen::BDCSVD<MatrixXd> svd(par_diff);
    const VectorXd s = svd.singularValues();

    const VectorXd s2 = s.cwiseProduct(s);
    const double sum_lambda  = s2.sum();
    const double sum_lambda2 = s2.cwiseProduct(s2).sum();

    double ess;
    if (sum_lambda2 <= 0.0 || sum_lambda <= 0.0) {
        // Degenerate (all-zero anomalies): treat as a single mode.
        ess = 1.0;
    } else {
        ess = (sum_lambda * sum_lambda) / sum_lambda2;
    }

    const double r_max = static_cast<double>(std::min(npar, nreal - 1));
    const double ess_ratio = (r_max > 0.0) ? (ess / r_max) : 1.0;

    return Rcpp::List::create(
        Rcpp::Named("ess")       = ess,
        Rcpp::Named("r_max")     = r_max,
        Rcpp::Named("ess_ratio") = ess_ratio
    );
}


//' Gaspari-Cohn Localisation Taper
//'
//' Evaluates the Gaspari & Cohn (1999) fifth-order piecewise-rational
//' compactly-supported correlation function on a matrix of distances. This
//' is the classical distance-based localisation taper: a smooth bump that
//' equals 1 at zero distance, decays to 0 at twice the localisation radius,
//' and is identically 0 beyond. Used to taper the Kalman gain when the
//' parameters and observations carry a spatial (or otherwise metric)
//' coordinate.
//'
//' With \eqn{z = d / c} (distance over localisation radius \eqn{c}):
//' \deqn{G(z) = \begin{cases}
//'   -\tfrac{1}{4}z^5 + \tfrac{1}{2}z^4 + \tfrac{5}{8}z^3 - \tfrac{5}{3}z^2 + 1
//'     & 0 \le z \le 1 \\
//'   \tfrac{1}{12}z^5 - \tfrac{1}{2}z^4 + \tfrac{5}{8}z^3 + \tfrac{5}{3}z^2
//'     - 5z + 4 - \tfrac{2}{3}z^{-1} & 1 < z \le 2 \\
//'   0 & z > 2.
//' \end{cases}}
//'
//' @param distances Matrix (npar x nobs). Non-negative parameter-to-
//'   observation distances.
//' @param radius Numeric scalar (> 0). Localisation radius \eqn{c}; the
//'   taper reaches 0 at distance \eqn{2c}.
//' @return Matrix (npar x nobs) of taper weights in \eqn{[0, 1]}.
//' @references
//' Gaspari, G. & Cohn, S.E. (1999). Construction of correlation functions
//' in two and three dimensions. \emph{Quarterly Journal of the Royal
//' Meteorological Society}, 125(554), 723--757.
//' @examples
//' d <- matrix(c(0, 0.5, 1, 1.5, 2, 3), nrow = 2L)
//' gaspari_cohn(d, radius = 1.0)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd gaspari_cohn(const Eigen::MatrixXd& distances, double radius)
{
    if (radius <= 0.0) {
        Rcpp::stop("`radius` must be a positive scalar.");
    }
    MatrixXd out(distances.rows(), distances.cols());
    for (int j = 0; j < distances.cols(); ++j) {
        for (int i = 0; i < distances.rows(); ++i) {
            const double d = distances(i, j);
            if (d < 0.0) {
                Rcpp::stop("`distances` must be non-negative.");
            }
            const double z = d / radius;
            double g;
            if (z <= 1.0) {
                g = -0.25 * std::pow(z, 5) + 0.5 * std::pow(z, 4)
                    + 0.625 * std::pow(z, 3) - (5.0 / 3.0) * z * z + 1.0;
            } else if (z <= 2.0) {
                g = (1.0 / 12.0) * std::pow(z, 5) - 0.5 * std::pow(z, 4)
                    + 0.625 * std::pow(z, 3) + (5.0 / 3.0) * z * z
                    - 5.0 * z + 4.0 - (2.0 / 3.0) / z;
            } else {
                g = 0.0;
            }
            // Guard the piecewise rounding at the knots.
            out(i, j) = std::min(1.0, std::max(0.0, g));
        }
    }
    return out;
}


//' Correlation-Based Automatic Localisation Taper
//'
//' Builds an \code{npar x nobs} localisation taper directly from the
//' ensemble, with no parameter / observation coordinates required. This is
//' the iterative-ensemble-smoother-native localisation of Luo & Bhakta
//' (2020): spurious sample correlations between a parameter and an
//' observation (an artefact of finite ensemble size) are damped, while
//' genuine correlations that stand above an estimated noise floor are
//' retained.
//'
//' The sample correlation \eqn{\rho_{ij}} between parameter-anomaly row
//' \eqn{i} and observation-anomaly row \eqn{j} is compared against a noise
//' floor \eqn{\theta}. When \code{threshold < 0} the floor is estimated by
//' destroying the parameter-observation link --- the realisation order of
//' the observation anomalies is randomly permuted, independently per
//' replicate, and the floor is taken as a high quantile (default 0.95) of
//' the resulting spurious \eqn{|\rho|} values. The permutation uses R's RNG,
//' so the estimate is reproducible under \code{set.seed()}.
//'
//' Two tapers are offered. \code{"hard"} keeps correlations above the floor
//' unchanged (weight 1) and zeroes the rest. \code{"soft"} applies a smooth,
//' monotone ramp \eqn{w_{ij} = \mathrm{clip}((|\rho_{ij}| - \theta) /
//' (1 - \theta), 0, 1)}, which downweights near-floor correlations rather
//' than thresholding them sharply.
//'
//' @param par_diff Matrix (npar x nreal). Parameter anomalies.
//' @param obs_diff Matrix (nobs x nreal). Observation anomalies.
//' @param threshold Numeric. Noise floor on \eqn{|\rho|}. Negative
//'   (default \code{-1}) triggers automatic estimation by permutation.
//' @param taper Character. \code{"hard"} (indicator) or \code{"soft"}
//'   (linear ramp above the floor).
//' @param n_shuffle Integer. Number of permutation replicates for the
//'   automatic floor (default 1; each replicate yields \code{npar * nobs}
//'   spurious samples). Ignored when \code{threshold >= 0}.
//' @param quantile Numeric in (0, 1). Quantile of the spurious-correlation
//'   distribution used as the floor (default 0.95). Ignored when
//'   \code{threshold >= 0}.
//' @return A list with \code{rho} (the npar x nobs taper), \code{threshold}
//'   (the floor used), \code{n_active} (count of entries with non-zero
//'   weight), and \code{frac_active} (that count over \code{npar * nobs}).
//' @references
//' Luo, X. & Bhakta, T. (2020). Automatic and adaptive localization for
//' ensemble-based history matching. \emph{Journal of Petroleum Science and
//' Engineering}, 184, 106559.
//' @examples
//' set.seed(1L)
//' npar <- 8L; nobs <- 5L; nreal <- 40L
//' pd <- matrix(rnorm(npar * nreal), npar, nreal)
//' # Make parameter 1 genuinely correlated with observation 1.
//' od <- matrix(rnorm(nobs * nreal), nobs, nreal)
//' od[1L, ] <- od[1L, ] + 2 * pd[1L, ]
//' loc <- correlation_localisation(pd, od)
//' loc$threshold
//' loc$rho[1L, 1L]
//' @export
// [[Rcpp::export]]
Rcpp::List correlation_localisation(
    const Eigen::MatrixXd& par_diff,
    const Eigen::MatrixXd& obs_diff,
    double threshold = -1.0,
    std::string taper = "hard",
    int n_shuffle = 1,
    double quantile = 0.95)
{
    const int npar  = par_diff.rows();
    const int nobs  = obs_diff.rows();
    const int nreal = par_diff.cols();
    if (obs_diff.cols() != nreal) {
        Rcpp::stop("`par_diff` and `obs_diff` must have the same number of "
                   "columns (realisations).");
    }
    if (nreal < 3) {
        Rcpp::stop("Correlation localisation needs at least 3 realisations.");
    }
    if (taper != "hard" && taper != "soft") {
        Rcpp::stop("`taper` must be \"hard\" or \"soft\".");
    }

    // Row norms (anomalies are assumed already mean-centred across columns).
    VectorXd par_norm(npar), obs_norm(nobs);
    for (int i = 0; i < npar; ++i) par_norm(i) = par_diff.row(i).norm();
    for (int j = 0; j < nobs; ++j) obs_norm(j) = obs_diff.row(j).norm();

    // Sample correlation matrix (npar x nobs).
    MatrixXd corr = par_diff * obs_diff.transpose();
    for (int j = 0; j < nobs; ++j) {
        for (int i = 0; i < npar; ++i) {
            const double denom = par_norm(i) * obs_norm(j);
            corr(i, j) = (denom > 0.0) ? (corr(i, j) / denom) : 0.0;
        }
    }

    // Estimate the noise floor by permuting the realisation order of the
    // observation anomalies, which destroys any true par-obs link and leaves
    // only finite-sample noise. R's RNG (set.seed-governed) drives the shuffle.
    double floor = threshold;
    if (threshold < 0.0) {
        if (n_shuffle < 1) n_shuffle = 1;
        std::vector<double> spurious;
        spurious.reserve(static_cast<std::size_t>(npar) * nobs * n_shuffle);
        std::vector<int> perm(nreal);
        for (int rep = 0; rep < n_shuffle; ++rep) {
            for (int k = 0; k < nreal; ++k) perm[k] = k;
            // Fisher-Yates using R's uniform generator.
            for (int k = nreal - 1; k > 0; --k) {
                const int m = static_cast<int>(unif_rand() * (k + 1));
                std::swap(perm[k], perm[m]);
            }
            MatrixXd obs_perm(nobs, nreal);
            for (int k = 0; k < nreal; ++k) obs_perm.col(k) = obs_diff.col(perm[k]);
            MatrixXd sc = par_diff * obs_perm.transpose();
            for (int j = 0; j < nobs; ++j) {
                for (int i = 0; i < npar; ++i) {
                    const double denom = par_norm(i) * obs_norm(j);
                    const double c = (denom > 0.0) ? (sc(i, j) / denom) : 0.0;
                    spurious.push_back(std::abs(c));
                }
            }
        }
        std::sort(spurious.begin(), spurious.end());
        // Linear-interpolation quantile on the sorted spurious magnitudes.
        const double q = std::min(std::max(quantile, 0.0), 1.0);
        const double pos = q * (static_cast<double>(spurious.size()) - 1.0);
        const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
        const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
        floor = spurious[lo] + (pos - lo) * (spurious[hi] - spurious[lo]);
    }

    // Build the taper.
    MatrixXd rho(npar, nobs);
    int n_active = 0;
    const double denom_soft = (1.0 - floor > 1e-12) ? (1.0 - floor) : 1e-12;
    for (int j = 0; j < nobs; ++j) {
        for (int i = 0; i < npar; ++i) {
            const double a = std::abs(corr(i, j));
            double w;
            if (taper == "hard") {
                w = (a >= floor) ? 1.0 : 0.0;
            } else {
                w = (a - floor) / denom_soft;
                w = std::min(1.0, std::max(0.0, w));
            }
            rho(i, j) = w;
            if (w > 0.0) ++n_active;
        }
    }

    const double frac_active =
        static_cast<double>(n_active) / (static_cast<double>(npar) * nobs);

    return Rcpp::List::create(
        Rcpp::Named("rho")         = rho,
        Rcpp::Named("threshold")   = floor,
        Rcpp::Named("n_active")    = n_active,
        Rcpp::Named("frac_active") = frac_active
    );
}


//' Localised Ensemble Solution Kernel (explicit-gain GLM form)
//'
//' Computes the IES Gauss-Levenberg-Marquardt update with state-space
//' covariance localisation applied as a Schur (elementwise) product on the
//' explicit Kalman gain. The standard SVD kernel [ensemble_solution()]
//' works in the reduced observation-anomaly subspace and never forms the
//' \code{npar x nobs} gain, so it cannot host localisation; this kernel
//' reconstructs the gain
//' \deqn{K = \Delta\Theta\, V\, \mathrm{diag}(s)\,
//'          \mathrm{diag}((s^2 + (\lambda+1))^{-1})\, U^{T},}
//' (with \eqn{U s V^{T}} the thin SVD of the weight-scaled observation
//' anomalies) tapers it as \eqn{K \circ \rho}, and applies it to the
//' weighted residuals.
//'
//' When \eqn{\rho \equiv 1} the result is identical (to truncation
//' tolerance) to [ensemble_solution()] with \code{use_approx = TRUE}; the
//' prior-scaling null-space correction (\code{upgrade_2}) is not part of the
//' localised path. The returned matrix follows the same sign convention as
//' [ensemble_solution()] --- it is the negative-direction step, applied to
//' the ensemble by subtraction (\code{par_new = par_old - upgrade}).
//'
//' @param par_diff Matrix (npar x nreal). Parameter anomalies.
//' @param obs_diff Matrix (nobs x nreal). Observation anomalies.
//' @param obs_resid Matrix (nobs x nreal). Observation residuals
//'   (sim - obs); see [ensemble_solution()] for the sign rationale.
//' @param weights Numeric vector (nobs). Observation weights
//'   (1 / sqrt(variance)).
//' @param rho Matrix (npar x nobs). Localisation taper in \eqn{[0, 1]},
//'   e.g. from [correlation_localisation()] or [gaspari_cohn()].
//' @param cur_lam Numeric. Current Marquardt lambda.
//' @param eigthresh Numeric. Eigenvalue truncation threshold (0-1).
//' @return Matrix (nreal x npar). Negative-direction parameter upgrade,
//'   applied by subtraction.
//' @references
//' Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
//' iterative ensemble smoother for efficient history matching and
//' uncertainty quantification. \emph{Computational Geosciences}, 17(4),
//' 689--703.
//' @examples
//' set.seed(1L)
//' npar <- 4L; nreal <- 20L; nobs <- 6L
//' par_diff  <- matrix(rnorm(npar * nreal), npar, nreal)
//' obs_diff  <- matrix(rnorm(nobs * nreal), nobs, nreal)
//' obs_resid <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
//' weights   <- rep(1, nobs)
//' rho       <- matrix(1, npar, nobs)          # no localisation
//' upg <- ensemble_solution_localised(
//'   par_diff, obs_diff, obs_resid, weights, rho, cur_lam = 1.0
//' )
//' dim(upg)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ensemble_solution_localised(
    Eigen::MatrixXd par_diff,
    Eigen::MatrixXd obs_diff,
    Eigen::MatrixXd obs_resid,
    Eigen::VectorXd weights,
    Eigen::MatrixXd rho,
    double cur_lam,
    double eigthresh = 1e-6)
{
    const int npar  = par_diff.rows();
    const int nobs  = obs_diff.rows();
    const int nreal = par_diff.cols();
    const double scale = 1.0 / std::sqrt(static_cast<double>(nreal - 1));

    if (rho.rows() != npar || rho.cols() != nobs) {
        Rcpp::stop("`rho` must be a npar x nobs matrix.");
    }

    DiagonalMatrix<double, Eigen::Dynamic> W = weights.asDiagonal();

    // Weight-scaled residuals and anomalies (mirrors ensemble_solution()).
    MatrixXd obs_resid_w = W * obs_resid;
    MatrixXd S = scale * (W * obs_diff);   // nobs x nreal
    MatrixXd P = scale * par_diff;          // npar x nreal

    // Thin SVD of the scaled observation anomalies.
    Eigen::BDCSVD<MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXd s = svd.singularValues();
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();

    int maxsing = s.size();
    if (eigthresh > 0) {
        const double thresh = eigthresh * s(0);
        for (int i = 0; i < s.size(); ++i) {
            if (s(i) < thresh) { maxsing = i; break; }
        }
        if (maxsing < s.size()) {
            s.conservativeResize(maxsing);
            U.conservativeResize(Eigen::NoChange, maxsing);
            V.conservativeResize(Eigen::NoChange, maxsing);
        }
    }

    VectorXd s2  = s.cwiseProduct(s);
    VectorXd ivec = (VectorXd::Ones(s2.size()) * (cur_lam + 1.0) + s2).cwiseInverse();

    // Explicit gain K = P V diag(s) diag(ivec) U^T  (npar x nobs).
    MatrixXd K = P * V * s.asDiagonal() * ivec.asDiagonal() * U.transpose();

    // Schur taper, then apply to the weighted residuals.
    K = K.cwiseProduct(rho);
    MatrixXd net = K * obs_resid_w;        // npar x nreal (net update direction)

    MatrixXd upgrade = -1.0 * net;         // negative-direction convention
    upgrade.transposeInPlace();            // nreal x npar
    return upgrade;
}
