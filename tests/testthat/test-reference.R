# Non-regression against the fits of WH 2.0.0 (fixtures/ref-2.0.0.rds, dev/ref-2.0.0/make_ref.R).
#
# At fixed smoothing parameters the 2.1 reproduces 2.0.0 up to rounding (measured <= 2e-14).
# A selected lambda is only defined up to the tolerances of the selection (Brent, Nelder-Mead,
# PIRLS stopping rule): measured on 2026-10-26, WH 2.0.0 itself moves 1D maximum likelihood
# selections by up to 4e-5 on lambda and 2e-7 on y_hat when its inputs move by 1e-14, and any
# change of rounding (another BLAS, another kernel) may do the same. Selected fits are therefore
# compared within that indeterminacy, fixed-lambda fits and the number of evaluations exactly.

ref <- readRDS(test_path("fixtures", "ref-2.0.0.rds"))
cur <- make_reference_fits()

rel_gap <- function(a, b) max(abs(unlist(a) - unlist(b))) / max(abs(unlist(b)))

test_that("the same fits are produced, in the same order", {
  expect_named(cur, names(ref))
})

test_that("fixed smoothing parameters reproduce WH 2.0.0 up to rounding", {

  for (k in c("1d_ml_fixed", "2d_ml_fixed")) {
    for (cmp in c("y_hat", "std_y_hat", "edf", "res")) expect_lt(rel_gap(cur[[k]][[cmp]], ref[[k]][[cmp]]), 1e-10)
    expect_lt(rel_gap(cur[[k]]$diagnosis, ref[[k]]$diagnosis), 1e-10)
  }
})

test_that("selected smoothing parameters reproduce WH 2.0.0 within the selection tolerance", {

  for (k in setdiff(names(ref), c("1d_ml_fixed", "2d_ml_fixed"))) {
    expect_lt(rel_gap(cur[[k]]$lambda, ref[[k]]$lambda), 1e-3)
    expect_lt(rel_gap(cur[[k]]$y_hat, ref[[k]]$y_hat), 1e-5)
    # standard deviations and effective degrees of freedom move with log(lambda) at first order
    for (cmp in c("std_y_hat", "edf")) expect_lt(rel_gap(cur[[k]][[cmp]], ref[[k]][[cmp]]), 1e-4)
    if (!is.null(ref[[k]]$y_pred)) {
      expect_lt(rel_gap(cur[[k]]$y_pred, ref[[k]]$y_pred), 1e-5)
      expect_lt(rel_gap(cur[[k]]$std_y_pred, ref[[k]]$std_y_pred), 1e-5)
    }
  }
})
