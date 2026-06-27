# Tests for the pestpp_available() capability probe. The probe must be a
# non-erroring boolean check so downstream code (examples, vignettes,
# conditional tests) can skip gracefully when no external binary exists.

test_that("pestpp_available() returns a single logical", {
  res <- pestpp_available()
  expect_type(res, "logical")
  expect_length(res, 1L)
  expect_false(is.na(res))
})

test_that("pestpp_available() never errors when a binary is absent", {
  # A name no PEST++ distribution ships -- must resolve to FALSE, not stop().
  expect_false(pestpp_available("pestpp-does-not-exist"))
})

test_that("pestpp_available() accepts the documented tool names", {
  for (nm in c("pestpp-ies", "pestpp-glm", "pestpp-swp", "pestpp-sen")) {
    res <- pestpp_available(nm)
    expect_type(res, "logical")
    expect_length(res, 1L)
  }
})

test_that("pestpp_available() rejects malformed input", {
  expect_error(pestpp_available(c("a", "b")))
  expect_error(pestpp_available(NA_character_))
  expect_error(pestpp_available(""))
  expect_error(pestpp_available(42L))
})
