# Non-regression against the fits of WH 2.0.0 (fixtures/ref-2.0.0.rds, dev/ref-2.0.0/make_ref.R).
#
# At fixed smoothing parameters the 2.1 reproduces 2.0.0 up to rounding (measured <= 2e-14).
#
# A selected lambda is only defined up to the tolerances of the selection (Brent, Nelder-Mead,
# PIRLS stopping rule), and a criterion flat around its optimum lets lambda move with no visible
# consequence on the fit. The selection is therefore judged by its impact on the total effective
# degrees of freedom sum(edf), in absolute terms whatever the size of the model (5 to 112 degrees
# of freedom across these fits): at most 1e-3 degree of freedom. Measured on 2026-10-26 by
# dev/ref-2.0.0/edf_gaps.R, largest gaps on 1d_ml_q1_BIC in both cases:
#   - WH 2.0.0 itself, inputs perturbed by 1e-14 (--calib): 2.3e-4 degree of freedom;
#   - the 2.1 against the reference: 2.6e-4 degree of freedom, the same order.
# The threshold leaves a factor of about 4 over both. In 2D the perturbation did not move the
# Nelder-Mead selections of 2.0.0 and the 2.1 reproduces them up to 3e-14 degree of freedom: the 2D
# margin is not calibrated, only inherited from 1D. Fitted values, which have no expression in
# degrees of freedom, keep a relative tolerance (largest gaps measured: 9.8e-7 on y_hat for the
# 2.1, 2e-7 for the perturbed 2.0.0). The number of evaluations, equal for all fits on 2026-10-26,
# is not checked: it follows the rounding through the comparisons of the optimiser.

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
    expect_lt(abs(sum(cur[[k]]$edf) - sum(ref[[k]]$edf)), 1e-3)
    expect_lt(rel_gap(cur[[k]]$y_hat, ref[[k]]$y_hat), 1e-5)
    if (!is.null(ref[[k]]$y_pred)) {
      expect_lt(rel_gap(cur[[k]]$y_pred, ref[[k]]$y_pred), 1e-5)
      expect_lt(rel_gap(cur[[k]]$std_y_pred, ref[[k]]$std_y_pred), 1e-5)
    }
  }
})
