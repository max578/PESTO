# Spectral Spread Effective Sample Size of a Parameter Ensemble

Diagnoses ensemble collapse (under-dispersion) by the participation
ratio of the parameter-anomaly covariance eigenspectrum. Given the
parameter anomalies \\\Delta\Theta\\ (deviations from the ensemble mean,
`npar x nreal`), the anomaly covariance is \\C = \Delta\Theta
\Delta\Theta^{T} / (N - 1)\\ with eigenvalues \\\lambda_i = s_i^2 / (N -
1)\\, where \\s_i\\ are the singular values of \\\Delta\Theta\\.

## Usage

``` r
ensemble_spread_ess(par_diff)
```

## Arguments

- par_diff:

  Matrix (npar x nreal). Parameter anomalies (deviations from the
  ensemble mean). At least 2 columns are required.

## Value

A list with components `ess` (the spectral spread-ESS), `r_max` (the
maximum attainable value \\\min(\mathrm{npar}, N - 1)\\), and
`ess_ratio` (`ess / r_max`, in \\(0, 1\]\\).

## Details

The spectral spread-ESS is the participation ratio \$\$\mathrm{ESS} =
\frac{(\sum_i \lambda_i)^2}{\sum_i \lambda_i^2} = \frac{(\sum_i
s_i^2)^2}{\sum_i s_i^4},\$\$ the effective number of directions carrying
variance. It is bounded in \\\[1, r\_{\max}\]\\ with \\r\_{\max} =
\min(\mathrm{npar}, N - 1)\\: equal to \\r\_{\max}\\ when variance is
spread isotropically across all modes, and approaching 1 when the
ensemble collapses onto a single direction. Because the ratio is
invariant to a global rescaling of the anomalies, it isolates the
*shape* of the collapse (directional degeneracy) from its *magnitude*;
magnitude is tracked separately by the R-side spread-retention ratio.

## References

Bretherton, C.S., Widmann, M., Dymnikov, V.P., Wallace, J.M. & Blade, I.
(1999). The effective number of spatial degrees of freedom of a
time-varying field. *Journal of Climate*, 12(7), 1990–2009.

## Examples

``` r
set.seed(1L)
# Healthy isotropic spread -> ESS near r_max
good <- matrix(rnorm(6L * 40L), 6L, 40L)
ensemble_spread_ess(good)$ess_ratio
#> [1] 0.8640052
# Collapsed onto one direction -> ESS near 1
v <- rnorm(6L)
bad <- outer(v, rnorm(40L)) + matrix(rnorm(6L * 40L, sd = 1e-3), 6L, 40L)
ensemble_spread_ess(bad)$ess_ratio
#> [1] 0.1666669
```
