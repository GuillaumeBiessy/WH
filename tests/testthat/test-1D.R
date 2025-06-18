data("portfolio_mort")
if (!interactive()) {
  pdf(NULL)
}

# Data----
d  <- portfolio_mort$d
ec <- portfolio_mort$ec
y <- log(d / ec)
y[d == 0] <- - 20
p <- length(d)

newdata <- 40:99

# Unit tests----

q <- 2
lambda <- 1e3
wt <- d
tXWz <- c(wt * y)

D <- diff(diag(p), differences = q)
P <- crossprod(D)
A <- lambda * P
diag(A) <- diag(A) + c(wt)

R <- chol(A)
K <- solve(R)
V <- chol2inv(R)

R_compact <- banded_to_compact(R, q)
A_compact <- banded_to_compact(A, q)

test_that("create P", {

  # Compact version
  expect_equal(
    banded_to_compact(P, q),
    create_P_compact_cpp(p, q),
    tolerance = .Machine$double.eps^.5
  )

  # Full version
  expect_equal(
    P, create_P_compact_cpp(p, q) |> compact_to_sym(),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("cholesky", {

  expect_equal(
    R,
    cholesky_compact_lapack(A_compact) |> compact_to_tri(),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("backsolve", {

  expect_equal(
    backsolve(R, backsolve(R, tXWz, transpose = TRUE)),
    backsolve_compact_cpp(R_compact, backsolve_compact_cpp(R_compact, tXWz, transpose = TRUE)),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("inverse cholesky", {

  expect_equal(
    K,
    invert_cholesky_compact_lapack(R_compact),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("chol2inv", {

  expect_equal(
    V,
    invert_cholesky_compact_lapack(R_compact) |> tcrossprod(),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("diag_V", {

  expect_equal(
    diag(V),
    diag_V_compact_cpp(R_compact),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("Regression models run and predict correctly", {

  for (q in 1:3) {
    expect_silent(WH_fixed_lambda(y = y, wt = wt, lambda = 1e2, q = q, verbose = 0))
    expect_silent(object <- WH_outer(y = y, wt = wt, q = q, verbose = 0))
    expect_silent(object_extra <- predict(object, newdata = newdata))
    expect_equal(
      sqrt(diag(vcov(object))),
      unname(object$std_y_hat),
      tolerance = .Machine$double.eps^.5)
    expect_equal(
      sqrt(diag(vcov(object_extra))),
      unname(object_extra$std_y_pred),
      tolerance = .Machine$double.eps^.5)
  }
})

test_that("Maximum Likelihood estimation runs and predicts correctly", {

  for (q in 1:3) {
    expect_silent(WH_fixed_lambda(d, ec, lambda = 1e2, q = q, verbose = 0))
    expect_silent(object <- WH_outer(d, ec, q = q, verbose = 0))
    expect_silent(object_extra <- predict(object, newdata = newdata))
    expect_equal(
      sqrt(diag(vcov(object))),
      unname(object$std_y_hat),
      tolerance = .Machine$double.eps^.5)
    expect_equal(
      sqrt(diag(vcov(object_extra))),
      unname(object_extra$std_y_pred),
      tolerance = .Machine$double.eps^.5
    )
  }
})
