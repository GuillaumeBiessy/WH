# Band kernels of the 2.1, against the equivalent dense algebra written in base R.

band_system <- function(d, lambda, q) {
  if (is.null(dim(d))) {
    P <- crossprod(diff(diag(length(d)), differences = q))
    band <- q
  } else {
    p <- dim(d)
    P <- lambda[[1]] / lambda[[2]] * (diag(p[[2]]) %x% crossprod(diff(diag(p[[1]]), differences = q[[1]]))) +
      crossprod(diff(diag(p[[2]]), differences = q[[2]])) %x% diag(p[[1]])
    lambda <- lambda[[2]]
    band <- q[[2]] * p[[1]]
  }
  A <- lambda * P
  diag(A) <- diag(A) + c(d)
  R <- chol(A)
  list(V = chol2inv(R), R_compact = banded_to_compact(R, band), band = band)
}

test_that("selected inverse: the band of A^{-1} from the band Cholesky factor (Takahashi)", {

  for (q in 1:3) {
    s <- band_system(portfolio_mort$d, 1e3, q)
    expect_equal(band_inverse_compact_cpp(s$R_compact), banded_to_compact(s$V, s$band), tolerance = 1e-10)
    expect_equal(diag_V_compact_cpp(s$R_compact), diag(s$V), tolerance = 1e-10)
  }
  s <- band_system(portfolio_LTC$d, c(1e1, 1e2), c(2, 2))
  expect_equal(band_inverse_compact_cpp(s$R_compact), banded_to_compact(s$V, s$band), tolerance = 1e-10)
  expect_equal(diag_V_compact_cpp(s$R_compact), diag(s$V), tolerance = 1e-10)
})

test_that("2D extrapolation variance: band C++ loop over new cells equals the dense formula", {

  fit <- WH(portfolio_LTC$d, portfolio_LTC$ec, verbose = 0)
  pred <- predict(fit, newdata = list(age = 65:104, duration = 0:17))
  # dense: V_22 = A V_11 A' + P_22^{-1}, A = - P_22^{-1} P_21, on the full grid
  n_tot <- c(40, 18)
  P <- fit$lambda[[1]] * (diag(n_tot[[2]]) %x% crossprod(diff(diag(n_tot[[1]]), differences = 2))) +
    fit$lambda[[2]] * (crossprod(diff(diag(n_tot[[2]]), differences = 2)) %x% diag(n_tot[[1]]))
  fit_pos <- pred$ind_fit
  P22_inv <- solve(P[- fit_pos, - fit_pos])
  A <- - P22_inv %*% P[- fit_pos, fit_pos]
  V22 <- A %*% vcov(fit) %*% t(A) + P22_inv
  expect_equal(c(pred$std_y_pred)[- fit_pos], sqrt(diag(V22)), tolerance = 1e-8)
})
