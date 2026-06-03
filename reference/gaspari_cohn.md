# Gaspari-Cohn Localisation Taper

Evaluates the Gaspari & Cohn (1999) fifth-order piecewise-rational
compactly-supported correlation function on a matrix of distances. This
is the classical distance-based localisation taper: a smooth bump that
equals 1 at zero distance, decays to 0 at twice the localisation radius,
and is identically 0 beyond. Used to taper the Kalman gain when the
parameters and observations carry a spatial (or otherwise metric)
coordinate.

## Usage

``` r
gaspari_cohn(distances, radius)
```

## Arguments

- distances:

  Matrix (npar x nobs). Non-negative parameter-to- observation
  distances.

- radius:

  Numeric scalar (\> 0). Localisation radius \\c\\; the taper reaches 0
  at distance \\2c\\.

## Value

Matrix (npar x nobs) of taper weights in \\\[0, 1\]\\.

## Details

With \\z = d / c\\ (distance over localisation radius \\c\\): \$\$G(z) =
\begin{cases} -\tfrac{1}{4}z^5 + \tfrac{1}{2}z^4 + \tfrac{5}{8}z^3 -
\tfrac{5}{3}z^2 + 1 & 0 \le z \le 1 \\ \tfrac{1}{12}z^5 -
\tfrac{1}{2}z^4 + \tfrac{5}{8}z^3 + \tfrac{5}{3}z^2 - 5z + 4 -
\tfrac{2}{3}z^{-1} & 1 \< z \le 2 \\ 0 & z \> 2. \end{cases}\$\$

## References

Gaspari, G. & Cohn, S.E. (1999). Construction of correlation functions
in two and three dimensions. *Quarterly Journal of the Royal
Meteorological Society*, 125(554), 723–757.

## Examples

``` r
d <- matrix(c(0, 0.5, 1, 1.5, 2, 3), nrow = 2L)
gaspari_cohn(d, radius = 1.0)
#>           [,1]       [,2] [,3]
#> [1,] 1.0000000 0.20833333    0
#> [2,] 0.6848958 0.01649306    0
```
