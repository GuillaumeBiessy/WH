# Convergence status published in the fit, and warnings emitted whatever verbose.

test_that("an ordinary selection converges, away from the boundary, without warning", {

  expect_no_warning(fit <- WH(portfolio_mort$d, portfolio_mort$ec, verbose = 0))
  expect_true(fit$convergence$converged)
  expect_false(fit$convergence$outer$boundary)
  expect_equal(fit$convergence$outer$evaluations, length(fit$lambda_track))
  expect_gt(fit$convergence$inner$iterations, 0)

  fit2 <- WH(portfolio_LTC$d, portfolio_LTC$ec, verbose = 0)
  expect_true(fit2$convergence$converged)
  expect_identical(fit2$convergence$outer$code, 0L)

  fixed <- WH(portfolio_mort$d, portfolio_mort$ec, lambda = 1e2, verbose = 0)
  expect_true(fixed$convergence$converged)
  expect_null(fixed$convergence$outer)
})

test_that("a smoothing parameter selected at the boundary is reported, even at verbose = 0", {

  # log-hazard exactly linear: lies in the null space of the penalty, the optimum is lambda = Inf
  x <- 50:94; ec <- rep(1e4, length(x)); d <- ec * exp(-9 + 0.08 * x); names(d) <- x
  expect_warning(fit <- WH(d, ec, verbose = 0), "upper boundary")
  expect_true(fit$convergence$outer$boundary)
  expect_equal(c(fit$y_hat), -9 + 0.08 * x, tolerance = 1e-4, ignore_attr = TRUE)

  # noise of variance 1e6 declared with variance 1 in the regression framework: REML wants to interpolate,
  # which takes lambda ~ 3e-9 at unit weights, below the search interval
  set.seed(1); y <- 1e3 * rnorm(length(x)); names(y) <- x
  expect_warning(fit <- WH(y = y, wt = rep(1, length(x)), verbose = 0), "lower boundary")
  expect_true(fit$convergence$outer$boundary)
})

test_that("a Nelder-Mead failure is reported and the fit published as not converged", {

  # a negative relative tolerance makes optim() stop on a degenerate simplex (code 10)
  expect_warning(fit <- WH(portfolio_LTC$d, portfolio_LTC$ec, tol_crit = - 1, verbose = 0), "did not converge")
  expect_false(fit$convergence$converged)
  expect_false(fit$convergence$outer$code == 0)
})

test_that("the PIRLS stops at its iteration cap (2.0.0 looped forever at tol_dev = 0)", {

  expect_warning(fit <- WH(portfolio_mort$d, portfolio_mort$ec, lambda = 1e2, tol_dev = 0, verbose = 0),
                 "PIRLS stopped after 200 iterations")
  expect_false(fit$convergence$converged)
  expect_equal(fit$convergence$inner$iterations, 200)
})
