// PESTO: Adaptive SVD and Linear Algebra Kernels
// Backends (CPU): randomised SVD, Apple Accelerate / LAPACK, Eigen.
//
// Planned future work: an optional GPU (cuSOLVER) backend. A flag-gated
// prototype was removed pending a complete implementation and is preserved in
// the project history.
//
// This module provides high-performance SVD computation with automatic
// backend selection based on matrix size and available hardware.
//
// References:
//   - Halko, Martinsson, Tropp (2011). Finding structure with randomness.
//   - Liberty et al. (2007). Randomized algorithms for the low-rank approximation.
//
// Copyright (c) 2026 Max Moldovan. Licensed under GPL-3 or any later version.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

// Platform-specific includes
// Note: We avoid including full Accelerate on macOS due to COMPLEX
// redefinition conflict with R headers. Instead, use LAPACK directly.
// The LAPACK routines are available via R's linked BLAS/LAPACK.
#define PESTO_HAS_ACCELERATE 0

// A GPU / cuSOLVER backend is planned future work; the prior flag-gated stub
// was removed (see the project history for the previous implementation).

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ============================================================================
// Randomised SVD (Halko-Martinsson-Tropp algorithm)
// State-of-the-art for large matrices where k << min(m,n)
// ============================================================================

//' Randomised SVD (Halko-Martinsson-Tropp Algorithm)
//'
//' Computes a rank-k approximation to the SVD using randomised
//' projections. This is asymptotically faster than full SVD for
//' problems where the target rank k is much smaller than min(m,n).
//'
//' Complexity: O(mn*k) vs O(mn*min(m,n)) for full SVD.
//'
//' @param A Matrix (m x n). Input matrix.
//' @param k Integer. Target rank (number of singular values to compute).
//' @param p Integer. Oversampling parameter (default 10). Higher = more accurate.
//' @param q Integer. Number of power iterations (default 2). Higher = better for
//'   matrices with slowly decaying singular values.
//' @return A list with components U (m x k), d (k), V (n x k).
//' @references
//' Halko, N., Martinsson, P.G., & Tropp, J.A. (2011). Finding structure
//' with randomness: Probabilistic algorithms for constructing approximate
//' matrix decompositions. SIAM Review, 53(2), 217-288.
//' @examples
//' set.seed(1L)
//' A <- matrix(rnorm(10 * 6), nrow = 10, ncol = 6)
//' res <- rsvd(A, k = 3L)
//' length(res$d)
//' A_hat <- res$u %*% diag(res$d) %*% t(res$v)
//' mean((A - A_hat)^2)
//' @export
// [[Rcpp::export]]
Rcpp::List rsvd(const Eigen::MatrixXd& A, int k, int p = 10, int q = 2) {
    int m = A.rows();
    int n = A.cols();
    int l = std::min(k + p, std::min(m, n));

    // Stage A: Form an approximate basis for the range of A
    // 1. Generate random Gaussian matrix
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    MatrixXd Omega(n, l);
    for (int j = 0; j < l; ++j)
        for (int i = 0; i < n; ++i)
            Omega(i, j) = dist(rng);

    // 2. Form Y = A * Omega
    MatrixXd Y = A * Omega;

    // 3. Power iterations for improved accuracy
    // (critical for matrices with slowly decaying singular values)
    for (int i = 0; i < q; ++i) {
        // QR factorisation for numerical stability
        Eigen::HouseholderQR<MatrixXd> qr1(Y);
        Y = qr1.householderQ() * MatrixXd::Identity(m, l);
        MatrixXd Z = A.transpose() * Y;
        Eigen::HouseholderQR<MatrixXd> qr2(Z);
        Z = qr2.householderQ() * MatrixXd::Identity(n, l);
        Y = A * Z;
    }

    // 4. QR of Y to get orthonormal basis Q
    Eigen::HouseholderQR<MatrixXd> qr(Y);
    MatrixXd Q = qr.householderQ() * MatrixXd::Identity(m, l);

    // Stage B: Form the small matrix B = Q^T * A and compute its SVD
    MatrixXd B = Q.transpose() * A;
    Eigen::BDCSVD<MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Extract rank-k approximation
    MatrixXd U = Q * svd.matrixU().leftCols(k);
    VectorXd d = svd.singularValues().head(k);
    MatrixXd V = svd.matrixV().leftCols(k);

    return Rcpp::List::create(
        Rcpp::Named("u") = U,
        Rcpp::Named("d") = d,
        Rcpp::Named("v") = V
    );
}

// ============================================================================
// Accelerate-backed SVD (macOS — uses Apple's vecLib/LAPACK via Accelerate)
// ============================================================================

// LAPACK SVD via R's linked BLAS/LAPACK
// Works on all platforms (macOS uses Accelerate/AMX under the hood)
extern "C" {
    void dgesvd_(const char* jobu, const char* jobvt,
                 const int* m, const int* n, double* a, const int* lda,
                 double* s, double* u, const int* ldu,
                 double* vt, const int* ldvt,
                 double* work, const int* lwork, int* info);
}

//' Hardware-Accelerated SVD via LAPACK
//'
//' Uses R's linked LAPACK (which on macOS is Apple Accelerate/AMX,
//' on Linux is typically OpenBLAS or MKL) for hardware-optimised
//' SVD computation.
//'
//' @param A Matrix (m x n). Input matrix.
//' @param thin Logical. If TRUE (default), compute thin SVD.
//' @return A list with components U, d, V.
//' @examples
//' set.seed(1L)
//' A <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
//' res <- accelerate_svd(A, thin = TRUE)
//' length(res$d)
//' all.equal(sort(res$d, decreasing = TRUE), svd(A)$d)
//' @export
// [[Rcpp::export]]
Rcpp::List accelerate_svd(const Eigen::MatrixXd& A, bool thin = true) {
    int m = A.rows();
    int n = A.cols();
    int k = std::min(m, n);

    // Copy to column-major for LAPACK
    std::vector<double> a_data(A.data(), A.data() + m * n);
    std::vector<double> s(k);
    std::vector<double> u(m * k);
    std::vector<double> vt(k * n);
    std::vector<double> work(1);
    int lwork = -1;
    int info = 0;

    char jobu = thin ? 'S' : 'A';
    char jobvt = thin ? 'S' : 'A';

    // Query optimal workspace size
    dgesvd_(&jobu, &jobvt, &m, &n, a_data.data(), &m,
            s.data(), u.data(), &m, vt.data(), &k,
            work.data(), &lwork, &info);

    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Compute SVD
    dgesvd_(&jobu, &jobvt, &m, &n, a_data.data(), &m,
            s.data(), u.data(), &m, vt.data(), &k,
            work.data(), &lwork, &info);

    if (info != 0) {
        Rcpp::stop("LAPACK dgesvd failed with info = %d", info);
    }

    // Map back to Eigen
    Eigen::Map<MatrixXd> U_mat(u.data(), m, k);
    Eigen::Map<VectorXd> s_vec(s.data(), k);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        Vt_mat(vt.data(), k, n);

    return Rcpp::List::create(
        Rcpp::Named("u") = MatrixXd(U_mat),
        Rcpp::Named("d") = VectorXd(s_vec),
        Rcpp::Named("v") = MatrixXd(Vt_mat.transpose())
    );
}


// ============================================================================
// Adaptive SVD: Automatically selects the best backend
// ============================================================================

//' Adaptive SVD with Automatic Backend Selection
//'
//' Automatically selects the SVD algorithm from the matrix dimensions:
//'
//' - **Randomised SVD** (Halko et al., 2011): when target rank k << min(m, n)
//' - **Apple Accelerate / LAPACK**: a full dense decomposition otherwise
//' - **Eigen BDCSVD**: cross-platform divide-and-conquer fallback
//'
//' @param A Matrix (m x n). Input matrix.
//' @param k Integer. Target rank. If 0 (default), computes full SVD.
//' @param method Character. Force a specific method: "auto" (default),
//'   "rsvd", "accelerate", "eigen".
//' @return A list with components U (m x k), d (k), V (n x k),
//'   plus `method_used` and `time_ms`.
//' @examples
//' set.seed(1L)
//' A <- matrix(rnorm(20 * 12), nrow = 20, ncol = 12)
//' res <- adaptive_svd(A, k = 5L, method = "auto")
//' length(res$d)
//' is.character(res$method_used)
//' @export
// [[Rcpp::export]]
Rcpp::List adaptive_svd(const Eigen::MatrixXd& A, int k = 0,
                         std::string method = "auto") {

    int m = A.rows();
    int n = A.cols();
    int mn = std::min(m, n);

    if (k <= 0 || k > mn) k = mn;

    auto t0 = std::chrono::high_resolution_clock::now();
    Rcpp::List result;
    std::string method_used;

    if (method == "auto") {
        // Decision logic:
        // 1. If k < 0.3 * min(m,n) and matrix is large: use randomised SVD
        // 2. On macOS: use Accelerate for full SVD
        // 3. Otherwise: Eigen BDCSVD

        if (k < static_cast<int>(0.3 * mn) && mn > 100) {
            method = "rsvd";
        } else {
            method = "accelerate"; // Uses R's linked LAPACK (Accelerate/MKL/OpenBLAS)
        }
    }

    if (method == "rsvd") {
        result = rsvd(A, k, 10, 2);
        method_used = "rsvd (Halko-Martinsson-Tropp)";
    } else if (method == "accelerate") {
        result = accelerate_svd(A, true);
        // Truncate to k
        MatrixXd U = Rcpp::as<MatrixXd>(result["u"]);
        VectorXd d = Rcpp::as<VectorXd>(result["d"]);
        MatrixXd V = Rcpp::as<MatrixXd>(result["v"]);
        result = Rcpp::List::create(
            Rcpp::Named("u") = U.leftCols(k),
            Rcpp::Named("d") = d.head(k),
            Rcpp::Named("v") = V.leftCols(k)
        );
        method_used = "LAPACK (platform-optimised)";
    } else if (method == "eigen") {
        // Eigen BDCSVD
        Eigen::BDCSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        result = Rcpp::List::create(
            Rcpp::Named("u") = MatrixXd(svd.matrixU().leftCols(k)),
            Rcpp::Named("d") = VectorXd(svd.singularValues().head(k)),
            Rcpp::Named("v") = MatrixXd(svd.matrixV().leftCols(k))
        );
        method_used = "Eigen BDCSVD";
    } else {
        Rcpp::stop("Unknown SVD method: %s. Use 'auto', 'rsvd', 'accelerate', or 'eigen'.",
                    method.c_str());
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    result["method_used"] = method_used;
    result["time_ms"] = ms;
    result["matrix_dims"] = Rcpp::IntegerVector::create(m, n);
    result["target_rank"] = k;

    return result;
}


// ============================================================================
// Ensemble solution using adaptive SVD backend selection
// ============================================================================

//' Ensemble Solution with Adaptive SVD Backend
//'
//' A variant of `ensemble_solution()` that selects the SVD backend
//' automatically -- a randomised SVD for low-rank problems, otherwise a dense
//' LAPACK / Accelerate decomposition -- and returns the upgrade together with
//' timing diagnostics. It is the convenient entry point when the best backend
//' for a given problem size is not known in advance. All computation is on the
//' CPU.
//'
//' @param par_diff Matrix (npar x nreal). Parameter anomalies.
//' @param obs_diff Matrix (nobs x nreal). Observation anomalies.
//' @param obs_resid Matrix (nobs x nreal). Observation residuals.
//' @param par_resid Matrix (npar x nreal). Parameter residuals.
//' @param weights Numeric vector (nobs). Observation weights.
//' @param parcov_inv Numeric vector (npar). Inverse parameter covariance diagonal.
//' @param Am Matrix. Random Am matrix for upgrade_2.
//' @param cur_lam Numeric. Marquardt lambda.
//' @param eigthresh Numeric. Eigenvalue truncation threshold.
//' @param use_approx Logical. Skip upgrade_2.
//' @param use_prior_scaling Logical. Scale by prior covariance.
//' @param iter Integer. Current iteration.
//' @param reg_factor Numeric. Regularisation factor.
//' @param svd_method Character. SVD method: "auto", "rsvd", "accelerate", "eigen".
//' @param target_rank Integer. Target rank for randomised SVD (0 = auto).
//' @return A list with upgrade matrix and performance diagnostics.
//' @examples
//' set.seed(1L)
//' npar  <- 4L
//' nreal <- 20L
//' nobs  <- 30L
//' par_diff  <- matrix(rnorm(npar * nreal), npar, nreal)
//' obs_diff  <- matrix(rnorm(nobs * nreal), nobs, nreal)
//' obs_resid <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
//' par_resid <- matrix(rnorm(npar * nreal, sd = 0.1), npar, nreal)
//' weights    <- rep(1, nobs)
//' parcov_inv <- rep(1, npar)
//' Am         <- matrix(0, 0, 0)
//' res <- ensemble_solution_adaptive(
//'   par_diff   = par_diff,
//'   obs_diff   = obs_diff,
//'   obs_resid  = obs_resid,
//'   par_resid  = par_resid,
//'   weights    = weights,
//'   parcov_inv = parcov_inv,
//'   Am         = Am,
//'   cur_lam    = 1.0,
//'   svd_method = "auto"
//' )
//' dim(res$upgrade)
//' @export
// [[Rcpp::export]]
Rcpp::List ensemble_solution_adaptive(
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
    double reg_factor = -1.0,
    std::string svd_method = "auto",
    int target_rank = 0)
{
    auto t_total = std::chrono::high_resolution_clock::now();

    int num_reals = par_diff.cols();
    double scale = 1.0 / std::sqrt(static_cast<double>(num_reals - 1));

    // Scale
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> W = weights.asDiagonal();
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> P = parcov_inv.asDiagonal();

    obs_resid = W * obs_resid;
    obs_diff = scale * (W * obs_diff);
    if (use_prior_scaling) {
        par_diff = scale * P * par_diff;
    } else {
        par_diff = scale * par_diff;
    }

    // Adaptive SVD on observation difference matrix
    if (target_rank <= 0) {
        target_rank = std::min(static_cast<int>(obs_diff.rows()),
                               static_cast<int>(obs_diff.cols()));
    }

    auto t_svd = std::chrono::high_resolution_clock::now();
    Rcpp::List svd_result = adaptive_svd(obs_diff, target_rank, svd_method);
    auto t_svd_end = std::chrono::high_resolution_clock::now();

    MatrixXd U = Rcpp::as<MatrixXd>(svd_result["u"]);
    VectorXd s = Rcpp::as<VectorXd>(svd_result["d"]);
    MatrixXd V = Rcpp::as<MatrixXd>(svd_result["v"]);
    std::string method_used = Rcpp::as<std::string>(svd_result["method_used"]);

    // Apply eigenvalue threshold
    int maxsing = s.size();
    if (eigthresh > 0 && s.size() > 0) {
        double thresh = eigthresh * s(0);
        for (int i = 0; i < s.size(); ++i) {
            if (s(i) < thresh) { maxsing = i; break; }
        }
        if (maxsing < s.size()) {
            s.conservativeResize(maxsing);
            U.conservativeResize(Eigen::NoChange, maxsing);
            V.conservativeResize(Eigen::NoChange, maxsing);
        }
    }

    // GLM update
    VectorXd s2 = s.cwiseProduct(s);
    VectorXd ivec = (VectorXd::Ones(s2.size()) * (cur_lam + 1.0) + s2).cwiseInverse();

    MatrixXd X1 = U.transpose() * obs_resid;
    MatrixXd X2 = ivec.asDiagonal() * X1;
    MatrixXd X3 = V * s.asDiagonal() * X2;
    MatrixXd upgrade_1 = -1.0 * par_diff * X3;
    upgrade_1.transposeInPlace();

    // Upgrade 2: null-space correction
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

    auto t_end = std::chrono::high_resolution_clock::now();
    double svd_ms = std::chrono::duration<double, std::milli>(t_svd_end - t_svd).count();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_total).count();

    return Rcpp::List::create(
        Rcpp::Named("upgrade") = upgrade_1,
        Rcpp::Named("svd_method") = method_used,
        Rcpp::Named("svd_time_ms") = svd_ms,
        Rcpp::Named("total_time_ms") = total_ms,
        Rcpp::Named("singular_values_used") = maxsing,
        Rcpp::Named("singular_values_total") = static_cast<int>(s.size())
    );
}
