# Inputs refused with an explicit message (WH 2.0.0 failed with cryptic errors, or not at all).

d <- portfolio_mort$d; ec <- portfolio_mort$ec

test_that("missing, infinite or negative observations are refused", {

  d_na <- d; d_na[3] <- NA
  expect_error(WH(d_na, ec, verbose = 0), "d should not contain missing or infinite values")
  ec_inf <- ec; ec_inf[3] <- Inf
  expect_error(WH(d, ec_inf, verbose = 0), "ec should not contain missing or infinite values")
  d_neg <- d; d_neg[3] <- - 1
  expect_error(WH(d_neg, ec, verbose = 0), "d should not contain negative values")
  expect_error(WH(y = log(d / ec), wt = - d, verbose = 0), "wt should not contain negative values")
})

test_that("events without exposure are refused", {

  ec_0 <- ec; ec_0[3] <- 0
  expect_error(WH(d, ec_0, verbose = 0), "ec should be positive wherever d is positive")
  expect_error(WH(d, ec_0, reg = TRUE, verbose = 0), "ec should be positive wherever d is positive")
})

test_that("a missing observation is refused on a positive weight, accepted on a zero weight", {

  y <- log(d / ec); y[3] <- NA
  expect_error(WH(y = y, wt = d, verbose = 0), "y should be finite wherever wt is positive")
  wt <- d; wt[3] <- 0
  expect_no_error(WH(y = y, wt = wt, verbose = 0))
})

test_that("smoothing parameters, penalization order and dimensions are checked", {

  expect_error(WH(d, ec, lambda = NA_real_, verbose = 0), "lambda should only have positive elements")
  expect_error(WH(d, ec, lambda = Inf, verbose = 0), "lambda should only have positive elements")
  expect_error(WH(d, ec, lambda = - 1, verbose = 0), "lambda should only have positive elements")
  expect_error(WH(d, ec, lambda = c(1, 1), verbose = 0), "lambda should be of length 1")
  expect_error(WH(d, ec, q = 1.5, verbose = 0), "q should be an integer vector")
  expect_error(WH(d[1:2], ec[1:2], q = 2, verbose = 0), "more than q data points")
  expect_error(WH(d, ec, criterion = "AICc", verbose = 0), "criterion should be one of")
  arr <- array(1, c(3, 3, 3))
  expect_error(WH(arr, arr, verbose = 0), "only handles smoothing in 1 or 2 dimensions")
})
