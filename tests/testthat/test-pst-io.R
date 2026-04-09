test_that("create_pest_scenario produces valid pesto_pst", {
  pars <- data.table::data.table(
    parnme = c("k1", "k2", "ss1"),
    partrans = "log",
    parchglim = "factor",
    parval1 = c(1.0, 0.5, 0.001),
    parlbnd = c(0.01, 0.001, 1e-6),
    parubnd = c(100, 50, 0.1),
    pargp = c("hydraulic", "hydraulic", "storage")
  )
  obs <- data.table::data.table(
    obsnme = c("h1", "h2", "h3", "h4"),
    obsval = c(10.5, 12.3, 8.7, 15.1),
    weight = c(1, 1, 0.5, 1),
    obgnme = c("heads", "heads", "heads", "heads")
  )

  pst <- create_pest_scenario(
    parameters = pars,
    observations = obs,
    model_command = "python run_model.py"
  )

  expect_s3_class(pst, "pesto_pst")
  expect_equal(pst$control_data$npar, 3L)
  expect_equal(pst$control_data$nobs, 4L)
  expect_equal(pst$control_data$npargp, 2L)
  expect_equal(pst$control_data$nobsgp, 1L)
})

test_that("write_pst and read_pst roundtrip", {
  pars <- data.table::data.table(
    parnme = c("k1", "k2"),
    partrans = c("log", "log"),
    parchglim = c("factor", "factor"),
    parval1 = c(1.0, 0.5),
    parlbnd = c(0.01, 0.001),
    parubnd = c(100, 50),
    pargp = c("hk", "hk"),
    scale = c(1.0, 1.0),
    offset = c(0.0, 0.0),
    dercom = c(1L, 1L)
  )
  obs <- data.table::data.table(
    obsnme = c("h1", "h2"),
    obsval = c(10.5, 12.3),
    weight = c(1.0, 1.0),
    obgnme = c("heads", "heads")
  )

  pst <- create_pest_scenario(pars, obs, "python model.py",
                               pestpp_options = list(
                                 ies_num_reals = 50,
                                 noptmax = 3
                               ))

  # Write and read back
  tmp <- tempfile(fileext = ".pst")
  on.exit(unlink(tmp))

  write_pst(pst, tmp)
  expect_true(file.exists(tmp))

  pst2 <- read_pst(tmp)
  expect_s3_class(pst2, "pesto_pst")
  expect_equal(pst2$control_data$npar, 2L)
  expect_equal(pst2$control_data$nobs, 2L)
})
