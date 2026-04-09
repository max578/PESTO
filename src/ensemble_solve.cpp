// PESTO: Scalable Tools for Estimation of Parameters
// Core ensemble solution kernel (ported from PEST++ EnsembleMethodUtils.cpp)
// Licensed under GPL-3

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::DiagonalMatrix;

//' Ensemble Solution Kernel (GLM form)
//'
//' Implements the core IES ensemble update equation using the
//' Gauss-Levenberg-Marquardt (GLM) formulation from Chen & Oliver (2013).
//' This is the computational hotspot of the iterative ensemble smoother.
//'
//' The update equation solves:
//' \deqn{\Delta\theta = -\Delta\theta' V s (s^2 + \lambda I)^{-1} U^T r}
//' where the SVD is performed on the scaled observation difference matrix.
//'
//' @param par_diff Matrix (npar x nreal). Parameter anomalies (deviations from mean).
//' @param obs_diff Matrix (nobs x nreal). Observation anomalies.
//' @param obs_resid Matrix (nobs x nreal). Observation residuals (obs - sim).
//' @param par_resid Matrix (npar x nreal). Parameter residuals (par - prior mean).
//' @param weights Numeric vector (nobs). Observation weights (1/sqrt(variance)).
//' @param parcov_inv Numeric vector (npar). Diagonal of inverse parameter covariance.
//' @param Am Matrix (npar x nreal-1). Random Am matrix for upgrade_2 (optional).
//' @param cur_lam Numeric. Current Marquardt lambda.
//' @param eigthresh Numeric. Eigenvalue truncation threshold (0-1).
//' @param use_approx Logical. If TRUE, skip upgrade_2 (prior-scaling correction).
//' @param use_prior_scaling Logical. Scale by prior covariance.
//' @param iter Integer. Current iteration number.
//' @param reg_factor Numeric. Regularisation factor for upgrade_2 blending.
//' @return Matrix (nreal x npar). Parameter upgrade vectors (one row per realisation).
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ensemble_solution(
    Eigen::MatrixXd par_diff,
    Eigen::MatrixXd obs_diff,
    Eigen::MatrixXd obs_resid,
    Eigen::MatrixXd par_resid,
    Eigen::VectorXd weights,
    Eigen::VectorXd parcov_inv,
    Eigen::MatrixXd Am,
    double cur_lam,
    double eigthresh = 1e-6,
    bool use_approx = true,
    bool use_prior_scaling = false,
    int iter = 1,
    double reg_factor = -1.0)
{
    int num_reals = par_diff.cols();
    double scale = 1.0 / std::sqrt(static_cast<double>(num_reals - 1));

    // Create diagonal matrices
    DiagonalMatrix<double, Eigen::Dynamic> W = weights.asDiagonal();
    DiagonalMatrix<double, Eigen::Dynamic> P = parcov_inv.asDiagonal();

    // Scale residuals by weights
    obs_resid = W * obs_resid;

    // Scale differences
    obs_diff = scale * (W * obs_diff);
    if (use_prior_scaling) {
        par_diff = scale * P * par_diff;
    } else {
        par_diff = scale * par_diff;
    }

    // SVD of scaled observation difference matrix
    // Using Eigen's BDCSVD for better performance on large matrices
    Eigen::BDCSVD<MatrixXd> svd(obs_diff, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXd s = svd.singularValues();
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();

    // Apply eigenvalue threshold
    int maxsing = s.size();
    if (eigthresh > 0) {
        double thresh = eigthresh * s(0);
        for (int i = 0; i < s.size(); ++i) {
            if (s(i) < thresh) {
                maxsing = i;
                break;
            }
        }
        if (maxsing < s.size()) {
            s.conservativeResize(maxsing);
            U.conservativeResize(Eigen::NoChange, maxsing);
            V.conservativeResize(Eigen::NoChange, maxsing);
        }
    }

    // Compute inverse term: (s^2 + (lambda+1)*I)^{-1}
    VectorXd s2 = s.cwiseProduct(s);
    VectorXd ivec = (VectorXd::Ones(s2.size()) * (cur_lam + 1.0) + s2).cwiseInverse();

    // Upgrade 1: standard GLM update
    MatrixXd X1 = U.transpose() * obs_resid;
    MatrixXd X2 = ivec.asDiagonal() * X1;
    MatrixXd X3 = V * s.asDiagonal() * X2;
    MatrixXd upgrade_1 = -1.0 * par_diff * X3;
    upgrade_1.transposeInPlace();

    // Upgrade 2: null-space correction (if not using approximation and iter > 1)
    if (!use_approx && iter > 1 && Am.rows() > 0 && Am.cols() > 0) {
        if (use_prior_scaling) {
            par_resid = P * par_resid;
        }
        MatrixXd x4 = Am.transpose() * par_resid;
        MatrixXd x5 = Am * x4;
        MatrixXd x6 = par_diff.transpose() * x5;
        MatrixXd x7 = V * ivec.asDiagonal() * V.transpose() * x6;
        MatrixXd upgrade_2;
        if (use_prior_scaling) {
            upgrade_2 = -1.0 * P * par_diff * x7;
        } else {
            upgrade_2 = -1.0 * (par_diff * x7);
        }
        if (reg_factor >= 0) {
            upgrade_1 = upgrade_1 + (reg_factor * upgrade_2.transpose());
        } else {
            upgrade_1 = upgrade_1 + upgrade_2.transpose();
        }
    }

    return upgrade_1;
}


//' Ensemble Solution Kernel (MDA / Evensen form)
//'
//' Implements the Multiple Data Assimilation (MDA) update from
//' Evensen (2018), Section 14.3.2. Uses low-rank representation
//' of the error covariance.
//'
//' @param par_diff Matrix (npar x nreal). Parameter anomalies.
//' @param obs_diff Matrix (nobs x nreal). Observation anomalies.
//' @param obs_resid Matrix (nobs x nreal). Observation residuals.
//' @param obs_err Matrix (nobs x nreal). Observation error realisations.
//' @param cur_lam Numeric. Inflation factor.
//' @param eigthresh Numeric. Eigenvalue truncation threshold.
//' @return Matrix (nreal x npar). Parameter upgrade vectors.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ensemble_solution_mda(
    Eigen::MatrixXd par_diff,
    Eigen::MatrixXd obs_diff,
    Eigen::MatrixXd obs_resid,
    Eigen::MatrixXd obs_err,
    double cur_lam = 1.0,
    double eigthresh = 1e-6)
{
    // Low rank Cee - Section 14.3.2 Evensen Book
    obs_err = obs_err.colwise() - obs_err.rowwise().mean();
    obs_err = std::sqrt(cur_lam) * obs_err;

    // First SVD on obs_diff
    Eigen::BDCSVD<MatrixXd> svd0(obs_diff, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXd s0 = svd0.singularValues();
    MatrixXd U0 = svd0.matrixU();

    // Apply threshold
    int maxsing = s0.size();
    if (eigthresh > 0) {
        double thresh = eigthresh * s0(0);
        for (int i = 0; i < s0.size(); ++i) {
            if (s0(i) < thresh) { maxsing = i; break; }
        }
        s0.conservativeResize(maxsing);
        U0.conservativeResize(Eigen::NoChange, maxsing);
    }

    MatrixXd s0_i = s0.cwiseInverse().asDiagonal();
    MatrixXd X0 = s0_i * U0.transpose() * obs_err;

    // Second SVD
    Eigen::BDCSVD<MatrixXd> svd1(X0, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXd s1 = svd1.singularValues();
    MatrixXd U1 = svd1.matrixU();

    VectorXd s1_2 = s1.cwiseProduct(s1);
    VectorXd s1_2i = (VectorXd::Ones(s1_2.size()) + s1_2).cwiseInverse();

    MatrixXd X1 = U0 * s0_i * U1;
    MatrixXd X4 = s1_2i.asDiagonal() * X1.transpose();
    MatrixXd X2 = X4 * obs_resid;
    MatrixXd X3 = X1 * X2;
    X3 = obs_diff.transpose() * X3;

    MatrixXd upgrade = -1.0 * par_diff * X3;
    upgrade.transposeInPlace();

    return upgrade;
}


//' Compute Phi (Objective Function) for Ensemble
//'
//' Calculates the weighted sum of squared residuals for each
//' realisation in the ensemble.
//'
//' @param residuals Matrix (nobs x nreal). Observation residuals.
//' @param weights Numeric vector (nobs). Observation weights.
//' @return Numeric vector (nreal). Phi value per realisation.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd compute_phi(
    const Eigen::MatrixXd& residuals,
    const Eigen::VectorXd& weights)
{
    int nreal = residuals.cols();
    VectorXd phi(nreal);

    VectorXd w2 = weights.cwiseProduct(weights);

    for (int i = 0; i < nreal; ++i) {
        VectorXd r = residuals.col(i);
        phi(i) = r.cwiseProduct(w2).dot(r);
    }

    return phi;
}
