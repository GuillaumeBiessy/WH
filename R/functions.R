# Misc----

build_D_mat <- function(p, q) {diff(diag(p), differences = q)}

compute_XZ_mat <- function(XZ) {map(XZ, \(X) {X[[2]] %x% X[[1]]}) |> reduce(cbind)}

compute_Xtheta <- function(X_mat, theta) {c(X_mat %*% (c(theta, recursive = TRUE)))}

compute_tXWz <- function(X_mat, z, wt) {c(t(X_mat) %*% c(wt * z))}

compute_tXWZ <- function(X_mat, Z_mat, wt) {t(X_mat) %*% (c(wt) * Z_mat)}

compute_diag_XPsitZ <- function(X_mat, Z_mat, PsiXZ) {rowSums((X_mat %*% PsiXZ) * Z_mat)}

cum_index <- function(n) {

  c(0, cumsum(n)[- length(n)]) |>
    map2(n, \(x,y) x + seq_len(y))
}

blockdiag <- function(L) {

  L  <- compact(L)
  n  <- length(L)

  n1 <- map(L, nrow)
  n2 <- map(L, ncol)

  i1 <- cum_index(n1)
  i2 <- cum_index(n2)

  M  <- matrix(0, do.call(sum, n1), do.call(sum, n2))
  for(i in seq_along(L)) M[i1[[i]], i2[[i]]] <- L[[i]]

  return(M)
}

compute_res_deviance <- function(D, D_hat) {

  D_diff <- D - D_hat
  log_D_diff <- ifelse(D == 0, 0, D * (log(D) - log(D_hat)))

  sign(D_diff) * sqrt(2 * (log_D_diff - D_diff))
}

compute_deviance <- function(D, D_hat) {

  res_deviance <- compute_res_deviance(D, D_hat)
  sum(res_deviance * res_deviance)
}

compute_inner_cond <- function(old_lambda, lambda) {

  max(abs(log(old_lambda) - log(lambda)))
}

compute_outer_cond <- function(D, old_D_hat, D_hat) {

  abs(compute_deviance(old_D_hat, D) / compute_deviance(D_hat, D) - 1)
}

#  Fit----

eigen_dec <- function(n, q) {

  D_mat <- build_D_mat(n, q)
  P  <- crossprod(D_mat)
  p  <- nrow(P)
  ei <- eigen(P, TRUE)
  U  <- ei$vectors[, order(ei$values)]

  X <- U[, seq_len(q), drop = FALSE]
  Z <- U[, q + seq_len(p - q), drop = FALSE]
  Sigma <- diag(sort(ei$values)[q + seq_len(p - q)])

  list(X = X, Z = Z, Sigma = Sigma)
}

build_var_comp_2d <- function(X, Z, Sigma) {

  q <- map(X, ncol)
  p <- map(Z, ncol)

  # Double interactions
  d12s <- diag(q[[2]]) %x% Sigma[[1]]
  d21s <- Sigma[[2]] %x% diag(q[[1]])

  d12d <- diag(p[[2]]) %x% Sigma[[1]]
  d21d <- Sigma[[2]] %x% diag(p[[1]])

  G_blocks <- list(d12s, d21s, d12d) |> map(~0 * .x) |> list() |> rep(2)
  G_blocks[[1]][c(1,3)] <- list(d12s, d12d)
  G_blocks[[2]][c(2,3)] <- list(d21s, d21d)

  dLambda <- G_blocks |> map_depth(2, diag) |> map(c, recursive = TRUE)

  build_G <- function(lambda) {

    Gm <- G_blocks |>
      map2(lambda, function(x, y) map(x, ~.x * y)) |>
      transpose() |>
      map(reduce, `+`) |>
      blockdiag()

    list(Gm = Gm, dG = 1 / diag(Gm))
  }

  list(dLambda = dLambda, build_G = build_G)
}

