# One block per bug of WH 2.0.0, each red on 2.0.0 and green once fixed.

test_that("regression framework: a zero-death cell does not turn the whole fit into NA", {

  d <- portfolio_mort$d; ec <- portfolio_mort$ec; d[5] <- 0
  fit <- WH(d, ec, reg = TRUE, verbose = 0)
  expect_false(anyNA(fit$y_hat))

  # the value of y on a zero weight is irrelevant: same fit as with any finite value there
  y <- log(d / ec); y[d == 0] <- - 20
  expect_equal(fit$y_hat, WH(y = y, wt = d, verbose = 0)$y_hat, tolerance = 1e-10)

  d2 <- portfolio_LTC$d; ec2 <- portfolio_LTC$ec
  fit2 <- WH(d2, ec2, reg = TRUE, lambda = c(1e1, 1e2), verbose = 0)
  expect_false(anyNA(fit2$y_hat))
})
