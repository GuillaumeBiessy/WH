data("portfolio_LTC")
if (!interactive()) {
  pdf(NULL)
}

# Data----
d  <- portfolio_LTC$d
ec <- portfolio_LTC$ec
y <- log(d / ec)
y[d == 0] <- - 20
p <- dim(d)

newdata <- list(age = 60:109, duration = 0:19)

# Unit tests----

q <- c(2, 2)
lambda <- c(1e1, 1e2)
wt <- d
tXWz <- c(wt * y)

D_eig <- list(diff(diag(p[[1]]), differences = q[[1]]), diff(diag(p[[2]]), differences = q[[2]]))
P <- lapply(D_eig, crossprod)
P_kro <- list(diag(p[[2]]) %x% P[[1]], P[[2]] %x% diag(p[[1]]))
sum_P_kro <- Reduce(`+`, P_kro)
lambda_P_kro <- mapply(`*`, lambda, P_kro, SIMPLIFY = FALSE)
A <- Reduce(`+`, lambda_P_kro)
diag(A) <- diag(A) + c(wt)

R <- chol(A)
K <- solve(R)
V <- chol2inv(R)

P_compact <- mapply(banded_to_compact, P, q, SIMPLIFY = FALSE)
R_compact <- banded_to_compact(R, q[[2]] * p[[1]])
A_compact <- banded_to_compact(A, q[[2]] * p[[1]])


test_that("create P", {

  # Compact version
  expect_equal(
    banded_to_compact(sum_P_kro, q[[2]] * p[[1]]),
    combine_P_compact(P_compact),
    tolerance = .Machine$double.eps^.5
  )

  # Full version
  expect_equal(
    sum_P_kro,
    combine_P_compact(P_compact) |> compact_to_sym(),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("P_y", {

  expect_equal(
    c(sum_P_kro %*% c(d)),
    get_prod_P_y_compact_cpp(d, P_compact),
    tolerance = .Machine$double.eps^.5
  )
})

test_that("P_K", {

  expect_equal(
    sum_P_kro %*% K,
    get_prod_P_K_compact_cpp(K, P_compact),
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
      c(object$std_y_hat),
      tolerance = .Machine$double.eps^.5)
    expect_equal(
      sqrt(diag(vcov(object_extra))),
      c(object_extra$std_y_pred),
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
      c(object$std_y_hat),
      tolerance = .Machine$double.eps^.5)
    expect_equal(
      sqrt(diag(vcov(object_extra))),
      c(object_extra$std_y_pred),
      tolerance = .Machine$double.eps^.5
    )
  }
})
