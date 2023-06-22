data("portfolio_mort")
if(!interactive()) pdf(NULL)

# Data----
d <- portfolio_mort$d
ec <- portfolio_mort$ec

y <- log(d / ec)
y[d == 0] <- - 20
wt <- d

compare_fits <- function(f1, f2, tolerance = if (edition_get() >= 3) testthat_tolerance()) {

  expect_equal(f1$y_hat, f2$y_hat, tolerance = tolerance)
  expect_equal(f1$std_y_hat, f2$std_y_hat, tolerance = 100 * tolerance)
  expect_equal(f1$diagnosis$REML, f2$diagnosis$REML, tolerance = tolerance)
}

compare_reml <- function(f1, f2, tolerance = if (edition_get() >= 3) testthat_tolerance()) {

  expect_equal(f1$diagnosis$REML, f2$diagnosis$REML, tolerance = tolerance)
}

# Regression----
test_that("Various way of invoking the regression framework are working", {

  ref_fixed_lambda <- WH_1d_fixed_lambda(y = y, wt = wt, lambda = 1e2, reg = TRUE)

  compare_fits(WH_1d(y = y, wt = wt, lambda = 1e2), ref_fixed_lambda)
  compare_fits(WH_1d(d, ec, framework = "reg", lambda = 1e2), ref_fixed_lambda)
  compare_fits(WH_1d(d, y = y, lambda = 1e2), ref_fixed_lambda)
})

test_that("Outer iteration method calls the right function and is the default method", {

  ref_outer <- WH_1d_outer(y = y, wt = wt, reg = TRUE)

  compare_fits(WH_1d(y = y, wt = wt, method = "outer"), ref_outer)
  compare_fits(WH_1d(y = y, wt = wt), ref_outer)
})

test_that("Performance iteration method calls the right function", {

  ref_perf <- WH_1d_perf(y = y, wt = wt, reg = TRUE)

  compare_fits(WH_1d(y = y, wt = wt, method = "perf"), ref_perf)
})

test_that("REML is the default criterion", {

  ref_outer <- WH_1d_outer(y = y, wt = wt, reg = TRUE)

  compare_fits(WH_1d(y = y, wt = wt, criterion = "REML"), ref_outer)
})

test_that("Outer and performance iteration methods give very close results", {

  ref_perf <- WH_1d_perf(y = y, wt = wt, reg = TRUE)
  ref_outer <- WH_1d_outer(y = y, wt = wt, reg = TRUE)

  compare_reml(ref_perf, ref_outer, tolerance = 1e-6)
})

test_that("Other smoothing parameter selection criteria are working as well", {

  compare_fits(WH_1d(y = y, wt = wt, criterion = "AIC"),
               WH_1d_outer(y = y, wt = wt, criterion = "AIC", reg = TRUE))
  compare_fits(WH_1d(y = y, wt = wt, criterion = "BIC"),
               WH_1d_outer(y = y, wt = wt, criterion = "BIC", reg = TRUE))
  compare_fits(WH_1d(y = y, wt = wt, criterion = "GCV"),
               WH_1d_outer(y = y, wt = wt, criterion = "GCV", reg = TRUE))
})

test_that("Rank reduction works", {

  ref_perf <- WH_1d_perf(y = y, wt = wt, reg = TRUE)
  ref_outer <- WH_1d_outer(y = y, wt = wt, reg = TRUE)
  ref_perf_red <- WH_1d_perf(y = y, wt = wt, p = 20, reg = TRUE)
  ref_outer_red <- WH_1d_outer(y = y, wt = wt, p = 20, reg = TRUE)

  compare_fits(WH_1d(y = y, wt = wt, method = "outer", p = 20), ref_outer_red)
  compare_fits(WH_1d(y = y, wt = wt, method = "perf", p = 20), ref_perf_red)
  compare_reml(ref_perf_red, ref_perf, tolerance = 1e-2)
  compare_reml(ref_outer_red, ref_outer, tolerance = 1e-2)
  compare_reml(ref_perf_red, ref_outer_red, tolerance = 1e-6)
})

# Maximum likelihood----

test_that("Supplying lambda calls the fixed lambda function directly in the ML framework", {
  compare_fits(WH_1d(d, ec, lambda = 1e2),
               WH_1d_fixed_lambda(d, ec, lambda = 1e2))
})

test_that("Outer iteration method calls the right function and is the default method", {

  ref_ml_outer <- WH_1d_outer(d, ec)

  compare_fits(WH_1d(d, ec, method = "outer"), ref_ml_outer)
  compare_fits(WH_1d(d, ec), ref_ml_outer)
})

test_that("Performance iteration method calls the right function", {

  ref_ml_perf <- WH_1d_perf(d, ec)

  compare_fits(WH_1d(d, ec, method = "perf"), ref_ml_perf)
})

test_that("REML is the default criterion", {

  ref_ml_outer <- WH_1d_outer(d, ec)

  compare_fits(WH_1d_outer(d, ec, criterion = "REML"), ref_ml_outer)
})

test_that("Outer and performance iteration methods give close enough results", {

  ref_ml_perf <- WH_1d_perf(d, ec)
  ref_ml_outer <- WH_1d_outer(d, ec)

  compare_fits(ref_ml_perf, ref_ml_outer, tolerance = 1e-2)
})

test_that("Other smoothing parameter selection criteria are working as well", {

  compare_fits(WH_1d(d, ec, criterion = "AIC"),
               WH_1d_outer(d, ec, criterion = "AIC"))
  compare_fits(WH_1d(d, ec, criterion = "BIC"),
               WH_1d_outer(d, ec, criterion = "BIC"))
  compare_fits(WH_1d(d, ec, criterion = "GCV"),
               WH_1d_outer(d, ec, criterion = "GCV"))
})

test_that("Rank reduction works", {

  ref_ml_perf <- WH_1d_perf(d, ec)
  ref_ml_outer <- WH_1d_outer(d, ec)
  ref_ml_perf_red <- WH_1d_perf(d, ec, p = 20)
  ref_ml_outer_red <- WH_1d_outer(d, ec, p = 20)

  compare_fits(WH_1d(d, ec, method = "perf", p = 20), ref_ml_perf_red)
  compare_fits(WH_1d(d, ec, method = "outer", p = 20), ref_ml_outer_red)
  compare_fits(ref_ml_perf_red, ref_ml_perf, tolerance = 1e-2)
  compare_fits(ref_ml_outer_red, ref_ml_outer, tolerance = 1e-2)
  compare_fits(ref_ml_perf_red, ref_ml_outer_red, tolerance = 1e-2)
})

# Plots----

test_that("Plot functions work", {

  expect_no_error({
    # Regression
    ref_perf <- WH_1d_perf(y = y, wt = wt, reg = TRUE)
    ref_outer <- WH_1d_outer(y = y, wt = wt, reg = TRUE)

    ref_perf |> plot()
    ref_outer |> plot()

    ref_perf |> plot("res")
    ref_outer |> plot("res")

    ref_perf |> plot("edf")
    ref_outer |> plot("edf")

    # Maximum likelihood
    ref_ml_perf <- WH_1d_perf(d, ec)
    ref_ml_outer <- WH_1d_outer(d, ec)

    ref_ml_perf |> plot()
    ref_ml_outer |> plot()

    ref_ml_perf |> plot("res")
    ref_ml_outer |> plot("res")

    ref_ml_perf |> plot("edf")
    ref_ml_outer |> plot("edf")
  })
})

# Extrapolation----

test_that("Extrapolation and extrapolation plots work", {

  newdata <- 18:99

  perf_extra_reg <- WH_1d_perf(y = y, wt = wt, reg = TRUE) |> predict(newdata)
  outer_extra_reg <- WH_1d_outer(y = y, wt = wt, reg = TRUE) |> predict(newdata)

  compare_fits(perf_extra_reg, outer_extra_reg, tolerance = 1e-6)

  perf_extra_ml <- WH_1d_perf(d, ec) |> predict(newdata)
  outer_extra_ml <- WH_1d_outer(d, ec) |> predict(newdata)

  compare_fits(perf_extra_ml, outer_extra_ml, tolerance = 1e-2)

  expect_no_error({
    perf_extra_reg |> plot()
    outer_extra_reg |> plot()

    perf_extra_ml |> plot()
    outer_extra_ml |> plot()
  })
})
