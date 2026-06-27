# Print and plot methods for PESTO ensemble-smoother results

[`print()`](https://rdrr.io/r/base/print.html) gives a one-screen
summary of a fitted iterative ensemble-smoother run (driver, problem
dimensions, the phi-convergence trace, the spread-ESS dispersion
diagnostic and the failure rate).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
objective-function (phi) convergence trace by delegating to
[`plot_phi()`](https://max578.github.io/PESTO/reference/plot_phi.md).
For posterior parameter distributions (prior vs posterior) use
[`plot_ensemble()`](https://max578.github.io/PESTO/reference/plot_ensemble.md).

## Usage

``` r
# S3 method for class 'pesto_ies_result'
print(x, ...)

# S3 method for class 'pesto_ies_result'
plot(x, ...)
```

## Arguments

- x:

  A `pesto_ies_result`, as returned by
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  or
  [`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md).

- ...:

  Further arguments. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) these are
  passed to
  [`plot_phi()`](https://max578.github.io/PESTO/reference/plot_phi.md)
  (e.g. `log_scale`, `title`);
  [`print()`](https://rdrr.io/r/base/print.html) ignores them.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot2` object.

## Examples

``` r
set.seed(1)
G <- matrix(rnorm(18L), 6L, 3L)
prior <- matrix(rnorm(150L), 50L, 3L, dimnames = list(NULL, paste0("p", 1:3)))
obs <- stats::setNames(as.numeric(G %*% c(1, -0.5, 2)) + rnorm(6L, sd = 0.05),
                       paste0("o", 1:6))
fit <- pesto_ies_callback(function(t) t %*% t(G), prior, obs,
                          obs_sd = 0.05, noptmax = 4L, verbose = FALSE)
print(fit)
#> $phi
#>      iteration realisation          phi
#>          <int>       <int>        <num>
#>   1:         1           1  7941.118063
#>   2:         1           2 28682.646430
#>   3:         1           3 21839.224359
#>   4:         1           4  1156.296969
#>   5:         1           5 17649.083100
#>  ---                                   
#> 196:         4          46     3.496086
#> 197:         4          47     3.519188
#> 198:         4          48     3.488991
#> 199:         4          49     3.540618
#> 200:         4          50     3.504251
#> 
#> $par_ensemble
#>     real_name        p1         p2       p3
#>        <char>     <num>      <num>    <num>
#>  1:    real_1 0.9977789 -0.4670468 2.016870
#>  2:    real_2 0.9983534 -0.4632015 2.016941
#>  3:    real_3 0.9977350 -0.4665664 2.016210
#>  4:    real_4 0.9975298 -0.4686381 2.017237
#>  5:    real_5 0.9964921 -0.4667302 2.016075
#>  6:    real_6 0.9920280 -0.4711418 2.014498
#>  7:    real_7 0.9963290 -0.4701987 2.015760
#>  8:    real_8 0.9964517 -0.4672419 2.016660
#>  9:    real_9 0.9955095 -0.4689989 2.015689
#> 10:   real_10 0.9935466 -0.4688744 2.015220
#> 11:   real_11 0.9950267 -0.4682778 2.015213
#> 12:   real_12 0.9962991 -0.4690119 2.015717
#> 13:   real_13 0.9980786 -0.4683072 2.016517
#> 14:   real_14 0.9955784 -0.4684583 2.015402
#> 15:   real_15 0.9976805 -0.4651971 2.017059
#> 16:   real_16 0.9944313 -0.4715011 2.014148
#> 17:   real_17 0.9942110 -0.4675184 2.015770
#> 18:   real_18 0.9949440 -0.4679142 2.014632
#> 19:   real_19 0.9959319 -0.4661065 2.015934
#> 20:   real_20 0.9955740 -0.4687642 2.015412
#> 21:   real_21 0.9979204 -0.4667084 2.016141
#> 22:   real_22 0.9975241 -0.4669806 2.016427
#> 23:   real_23 0.9946814 -0.4696499 2.014124
#> 24:   real_24 0.9968656 -0.4653692 2.017281
#> 25:   real_25 0.9973013 -0.4655737 2.015390
#> 26:   real_26 0.9972701 -0.4663145 2.016145
#> 27:   real_27 0.9954160 -0.4654056 2.015303
#> 28:   real_28 0.9949028 -0.4674397 2.015218
#> 29:   real_29 0.9967721 -0.4698751 2.017435
#> 30:   real_30 0.9970485 -0.4686842 2.016182
#> 31:   real_31 0.9946093 -0.4708673 2.014420
#> 32:   real_32 0.9966248 -0.4688093 2.014903
#> 33:   real_33 0.9965627 -0.4689012 2.016337
#> 34:   real_34 0.9950489 -0.4682658 2.015683
#> 35:   real_35 0.9959731 -0.4697165 2.015565
#> 36:   real_36 0.9938625 -0.4685599 2.014711
#> 37:   real_37 0.9975207 -0.4688102 2.015241
#> 38:   real_38 1.0001103 -0.4634095 2.016747
#> 39:   real_39 0.9962975 -0.4664885 2.016898
#> 40:   real_40 0.9945953 -0.4668907 2.015285
#> 41:   real_41 0.9967224 -0.4671763 2.015273
#> 42:   real_42 0.9976418 -0.4641578 2.018089
#> 43:   real_43 0.9999711 -0.4677257 2.017319
#> 44:   real_44 0.9956302 -0.4690057 2.015605
#> 45:   real_45 0.9985686 -0.4643660 2.017741
#> 46:   real_46 0.9960891 -0.4690827 2.016501
#> 47:   real_47 0.9944259 -0.4690015 2.015028
#> 48:   real_48 0.9970617 -0.4681392 2.017769
#> 49:   real_49 0.9926891 -0.4697880 2.014754
#> 50:   real_50 0.9978322 -0.4680048 2.015447
#>     real_name        p1         p2       p3
#>        <char>     <num>      <num>    <num>
#> 
#> $obs_ensemble
#>     real_name        o1        o2       o3       o4         o5        o6
#>        <char>     <num>     <num>    <num>    <num>      <num>     <num>
#>  1:    real_1 -2.105676 -4.628359 1.166150 1.643743 -0.4099503 0.9028742
#>  2:    real_2 -2.104206 -4.625571 1.167964 1.643482 -0.4039489 0.9039687
#>  3:    real_3 -2.105004 -4.626550 1.165720 1.643556 -0.4092279 0.9024739
#>  4:    real_4 -2.106523 -4.630391 1.165854 1.643815 -0.4124440 0.9028039
#>  5:    real_5 -2.104222 -4.626601 1.166514 1.641629 -0.4098829 0.9033031
#>  6:    real_6 -2.102596 -4.627185 1.165930 1.635926 -0.4179976 0.9037573
#>  7:    real_7 -2.105615 -4.628494 1.164298 1.642442 -0.4151751 0.9017874
#>  8:    real_8 -2.104809 -4.628282 1.166911 1.641695 -0.4106792 0.9036887
#>  9:    real_9 -2.104472 -4.627602 1.165594 1.640772 -0.4136302 0.9028604
#> 10:   real_10 -2.102891 -4.626831 1.166778 1.637624 -0.4140811 0.9040766
#> 11:   real_11 -2.103523 -4.626104 1.165878 1.639803 -0.4126914 0.9030886
#> 12:   real_12 -2.104991 -4.627527 1.164958 1.642034 -0.4133900 0.9022336
#> 13:   real_13 -2.106259 -4.628451 1.164776 1.644622 -0.4117512 0.9018030
#> 14:   real_14 -2.104074 -4.626554 1.165525 1.640730 -0.4127855 0.9027438
#> 15:   real_15 -2.104830 -4.627430 1.167510 1.643013 -0.4071895 0.9038541
#> 16:   real_16 -2.104059 -4.626233 1.163320 1.639885 -0.4177433 0.9013147
#> 17:   real_17 -2.102987 -4.626925 1.167622 1.638245 -0.4118211 0.9045790
#> 18:   real_18 -2.102933 -4.624564 1.165502 1.639586 -0.4121595 0.9027498
#> 19:   real_19 -2.103479 -4.625930 1.167182 1.640551 -0.4091222 0.9038723
#> 20:   real_20 -2.104226 -4.626803 1.165364 1.640815 -0.4132496 0.9026377
#> 21:   real_21 -2.105147 -4.626468 1.165407 1.643898 -0.4093803 0.9022018
#> 22:   real_22 -2.105209 -4.627374 1.165902 1.643336 -0.4099271 0.9026901
#> 23:   real_23 -2.103299 -4.624768 1.164151 1.639720 -0.4148619 0.9018090
#> 24:   real_24 -2.104542 -4.628199 1.168342 1.641755 -0.4077217 0.9046655
#> 25:   real_25 -2.103739 -4.624080 1.165732 1.642598 -0.4078567 0.9024428
#> 26:   real_26 -2.104550 -4.626305 1.166181 1.642740 -0.4089991 0.9028924
#> 27:   real_27 -2.102422 -4.624110 1.167307 1.639543 -0.4082225 0.9039734
#> 28:   real_28 -2.103040 -4.625519 1.166469 1.639349 -0.4114653 0.9035216
#> 29:   real_29 -2.106775 -4.631882 1.165998 1.642975 -0.4145669 0.9031301
#> 30:   real_30 -2.105589 -4.628179 1.165044 1.643109 -0.4126552 0.9021857
#> 31:   real_31 -2.104031 -4.626336 1.163843 1.639963 -0.4167309 0.9016729
#> 32:   real_32 -2.104590 -4.625515 1.163887 1.642528 -0.4129632 0.9012770
#> 33:   real_33 -2.105487 -4.628770 1.165499 1.642393 -0.4131458 0.9026456
#> 34:   real_34 -2.103823 -4.627132 1.166394 1.639813 -0.4126735 0.9035185
#> 35:   real_35 -2.105036 -4.627772 1.164654 1.641736 -0.4145602 0.9020834
#> 36:   real_36 -2.102619 -4.625414 1.166123 1.638054 -0.4134932 0.9034597
#> 37:   real_37 -2.105362 -4.626100 1.163518 1.643943 -0.4126748 0.9008608
#> 38:   real_38 -2.105287 -4.624972 1.166158 1.646357 -0.4036813 0.9022629
#> 39:   real_39 -2.104493 -4.628279 1.167740 1.641208 -0.4095949 0.9043329
#> 40:   real_40 -2.102621 -4.625318 1.167118 1.638688 -0.4107376 0.9040511
#> 41:   real_41 -2.104085 -4.625111 1.165162 1.642169 -0.4104684 0.9021829
#> 42:   real_42 -2.104939 -4.628949 1.169299 1.642587 -0.4056477 0.9052627
#> 43:   real_43 -2.107659 -4.629452 1.164433 1.647427 -0.4102616 0.9012345
#> 44:   real_44 -2.104499 -4.627397 1.165394 1.640970 -0.4135992 0.9026788
#> 45:   real_45 -2.105405 -4.628163 1.168013 1.644145 -0.4056514 0.9040929
#> 46:   real_46 -2.105380 -4.629354 1.165974 1.641686 -0.4135790 0.9031179
#> 47:   real_47 -2.103384 -4.626339 1.165755 1.639074 -0.4139805 0.9031246
#> 48:   real_48 -2.106318 -4.631288 1.167132 1.642892 -0.4118527 0.9038849
#> 49:   real_49 -2.102509 -4.626631 1.166444 1.636556 -0.4157372 0.9039840
#> 50:   real_50 -2.105293 -4.625905 1.163953 1.644184 -0.4113580 0.9011139
#>     real_name        o1        o2       o3       o4         o5        o6
#>        <char>     <num>     <num>    <num>    <num>      <num>     <num>
#> 
#> $iterations
#> $iterations[[1]]
#> $iterations[[1]]$lambda
#> [1] 1
#> 
#> $iterations[[1]]$mean_phi
#> [1] 20268.07
#> 
#> $iterations[[1]]$n_failures
#> [1] 0
#> 
#> $iterations[[1]]$spread_ess
#> [1] 1.600024
#> 
#> $iterations[[1]]$spread_ess_ratio
#> [1] 0.5333414
#> 
#> $iterations[[1]]$inflation_method
#> [1] "none"
#> 
#> $iterations[[1]]$inflation_factor
#> [1] 1
#> 
#> $iterations[[1]]$retention
#> [1] NA
#> 
#> $iterations[[1]]$localisation
#> [1] "none"
#> 
#> $iterations[[1]]$loc_threshold
#> [1] NA
#> 
#> $iterations[[1]]$loc_frac_active
#> [1] NA
#> 
#> 
#> $iterations[[2]]
#> $iterations[[2]]$lambda
#> [1] 1
#> 
#> $iterations[[2]]$mean_phi
#> [1] 3.510204
#> 
#> $iterations[[2]]$n_failures
#> [1] 0
#> 
#> $iterations[[2]]$spread_ess
#> [1] 1.601765
#> 
#> $iterations[[2]]$spread_ess_ratio
#> [1] 0.5339218
#> 
#> $iterations[[2]]$inflation_method
#> [1] "none"
#> 
#> $iterations[[2]]$inflation_factor
#> [1] 1
#> 
#> $iterations[[2]]$retention
#> [1] NA
#> 
#> $iterations[[2]]$localisation
#> [1] "none"
#> 
#> $iterations[[2]]$loc_threshold
#> [1] NA
#> 
#> $iterations[[2]]$loc_frac_active
#> [1] NA
#> 
#> 
#> $iterations[[3]]
#> $iterations[[3]]$lambda
#> [1] 1
#> 
#> $iterations[[3]]$mean_phi
#> [1] 3.510128
#> 
#> $iterations[[3]]$n_failures
#> [1] 0
#> 
#> $iterations[[3]]$spread_ess
#> [1] 1.603498
#> 
#> $iterations[[3]]$spread_ess_ratio
#> [1] 0.5344993
#> 
#> $iterations[[3]]$inflation_method
#> [1] "none"
#> 
#> $iterations[[3]]$inflation_factor
#> [1] 1
#> 
#> $iterations[[3]]$retention
#> [1] NA
#> 
#> $iterations[[3]]$localisation
#> [1] "none"
#> 
#> $iterations[[3]]$loc_threshold
#> [1] NA
#> 
#> $iterations[[3]]$loc_frac_active
#> [1] NA
#> 
#> 
#> $iterations[[4]]
#> $iterations[[4]]$lambda
#> [1] 1
#> 
#> $iterations[[4]]$mean_phi
#> [1] 3.510053
#> 
#> $iterations[[4]]$n_failures
#> [1] 0
#> 
#> $iterations[[4]]$spread_ess
#> [1] 1.605222
#> 
#> $iterations[[4]]$spread_ess_ratio
#> [1] 0.535074
#> 
#> $iterations[[4]]$inflation_method
#> [1] "none"
#> 
#> $iterations[[4]]$inflation_factor
#> [1] 1
#> 
#> $iterations[[4]]$retention
#> [1] NA
#> 
#> $iterations[[4]]$localisation
#> [1] "none"
#> 
#> $iterations[[4]]$loc_threshold
#> [1] NA
#> 
#> $iterations[[4]]$loc_frac_active
#> [1] NA
#> 
#> 
#> 
#> $runtime_seconds
#> [1] 0.002
#> 
#> $n_forward_evals
#> [1] 250
#> 
#> $failure_rate
#> [1] 0
#> 
#> $converged
#> [1] FALSE
#> 
#> $n_iterations
#> [1] 4
#> 
#> $fidelity
#> NULL
#> 
#> $obs_target
#>         o1         o2         o3         o4         o5         o6 
#> -2.1198695 -4.6045419  1.2417415  1.6633979 -0.4359134  0.8684248 
#> 
#> $obs_sd
#>   o1   o2   o3   o4   o5   o6 
#> 0.05 0.05 0.05 0.05 0.05 0.05 
#> 
#> $weights
#> o1 o2 o3 o4 o5 o6 
#> 20 20 20 20 20 20 
#> 
#> attr(,"class")
#> [1] "pesto_ies_callback_result" "pesto_ies_result"         
p <- plot(fit)        # phi convergence (a ggplot2 object)
```
