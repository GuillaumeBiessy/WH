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

test_that("the starting smoothing parameters are found whatever the scale of the weights", {

  wt <- c(portfolio_LTC$d)
  P <- get_model_matrices(dim(portfolio_LTC$d), c(2, 2))$P
  n <- dim(portfolio_LTC$d)
  diag_P <- list(rep(P[[1]][nrow(P[[1]]), ], n[[2]]), rep(P[[2]][nrow(P[[2]]), ], each = n[[1]]))
  ref <- find_starting_lambda(wt, diag_P, 0.2)
  # the root for c * wt is c times the root for wt; 2.0.0 failed outside [1e-6, 1e12]
  for (scale in c(1e-12, 1e-6, 1e6, 1e12)) {
    expect_equal(find_starting_lambda(wt * scale, diag_P, 0.2) / scale, ref, tolerance = 1e-4)
  }
  y <- log(portfolio_LTC$d / portfolio_LTC$ec)
  expect_error(WH(y = y, wt = portfolio_LTC$d, lambda_start_ratio = 1, verbose = 0), "strictly between 0 and 1")
  expect_error(WH(y = y, wt = portfolio_LTC$d, lambda_start_ratio = 0, verbose = 0), "strictly between 0 and 1")
})

test_that("an evaluation whose Cholesky factorization fails does not stop the selection", {

  # weights x 1e-8: Nelder-Mead explores lambda where W + lambda P is numerically singular (2.0.0: error)
  y <- log(portfolio_LTC$d / portfolio_LTC$ec)
  expect_no_error(fit <- WH(y = y, wt = portfolio_LTC$d * 1e-8, verbose = 0))
  expect_true(all(is.finite(fit$lambda)) && is.finite(fit$diagnosis$REML))
})

test_that("predict adds the observed positions to newdata and handles non contiguous positions", {

  fit <- WH(portfolio_mort$d, portfolio_mort$ec, verbose = 0)
  ref <- predict(fit, newdata = 40:110)

  # newdata beyond the observations only: the observed positions are added (2.0.0 recycled silently)
  pred <- predict(fit, newdata = 100:110)
  expect_equal(names(pred$y_pred), as.character(c(50:94, 100:110)))
  expect_equal(pred$y_pred, ref$y_pred[names(pred$y_pred)])
  expect_equal(unname(pred$y_pred[as.character(50:94)]), unname(fit$y_hat), tolerance = 1e-6)

  # non contiguous newdata: the value at 105 is the value at 105 (2.0.0 returned the value at 100)
  pred <- predict(fit, newdata = c(45:49, 105))
  expect_equal(names(pred$y_pred), as.character(c(45:94, 105)))
  expect_equal(pred$y_pred, ref$y_pred[names(pred$y_pred)])
  expect_equal(pred$std_y_pred, ref$std_y_pred[names(pred$y_pred)])
  expect_equal(unname(sqrt(diag(vcov(pred)))), unname(pred$std_y_pred), tolerance = 1e-8)

  expect_error(predict(fit, newdata = 100.5), "grid of the observations")

  fit2 <- WH(portfolio_LTC$d, portfolio_LTC$ec, verbose = 0)
  ref2 <- predict(fit2, newdata = list(age = 70:109, duration = 0:19))  # the grid spanned by pred2
  pred2 <- predict(fit2, newdata = list(age = c(105, 109), duration = 19))
  rows <- rownames(pred2$y_pred); cols <- colnames(pred2$y_pred)
  expect_equal(rows, as.character(c(70:99, 105, 109)))
  expect_equal(pred2$y_pred[as.character(70:99), as.character(0:14)], fit2$y_hat)
  expect_equal(cols, as.character(c(0:14, 19)))
  expect_equal(pred2$y_pred, ref2$y_pred[rows, cols])
  expect_equal(pred2$std_y_pred, ref2$std_y_pred[rows, cols])
  expect_equal(sqrt(diag(vcov(pred2))), c(pred2$std_y_pred), tolerance = 1e-8)
})

test_that("the 1D log pseudo-determinant of the penalty is exact, also for long series", {

  # exact integer determinant of D D' by fraction-free (Bareiss) elimination
  bareiss <- function(M) {
    prev <- 1
    for (k in 1:(nrow(M) - 1)) {
      for (i in (k + 1):nrow(M)) for (j in (k + 1):nrow(M)) M[i, j] <- (M[i, j] * M[k, k] - M[i, k] * M[k, j]) / prev
      prev <- M[k, k]
    }
    M[nrow(M), nrow(M)]
  }
  for (q in 1:6) for (n in (q + 2):(q + 9)) {
    G <- tcrossprod(diff(diag(n), differences = q))
    expect_equal(log_pdet_P(n, q), log(bareiss(G)), tolerance = 1e-12)
  }
  # 2.0.0 summed the logs of LAPACK eigenvalues: NaN here (smallest nonzero eigenvalue computed negative)
  set.seed(1); d <- rpois(3000, 50); names(d) <- seq_along(d)
  fit <- WH(d, rep(1000, 3000), q = 3, lambda = 1e6, verbose = 0)
  expect_true(is.finite(fit$diagnosis$REML))
})
