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

test_that("data without names get default positions, in 1D and in 2D", {

  d <- unname(portfolio_mort$d); ec <- unname(portfolio_mort$ec)
  expect_warning(capture.output(fit <- WH(d, ec)), "No names found for d")
  expect_equal(names(fit$y_hat), as.character(seq_along(d)))
  expect_output(print(fit), paste0("1  to  ", length(d)))
  expect_no_warning(pred <- predict(fit, newdata = seq_len(length(d) + 5)))
  expect_length(pred$y_pred, length(d) + 5)

  expect_warning(capture.output(fit_reg <- WH(y = log(d / ec), wt = d)), "No names found for y")
  expect_equal(names(fit_reg$y_hat), as.character(seq_along(d)))

  d2 <- unname(portfolio_LTC$d); ec2 <- unname(portfolio_LTC$ec)
  fit2 <- suppressWarnings(WH(d2, ec2, lambda = c(1e1, 1e2), verbose = 0))
  expect_equal(dimnames(fit2$y_hat), list(as.character(seq_len(nrow(d2))), as.character(seq_len(ncol(d2)))))
  expect_equal(unname(fit2$y_hat), unname(WH(portfolio_LTC$d, portfolio_LTC$ec, lambda = c(1e1, 1e2), verbose = 0)$y_hat))
})
