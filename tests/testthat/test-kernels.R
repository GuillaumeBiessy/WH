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
