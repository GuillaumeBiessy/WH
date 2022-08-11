build_D_mat <- function(n, q) {diff(diag(n), differences = q)}

compute_XZ_mat <- function(XZ) {
  lapply(XZ, \(X) {X[[2]] %x% X[[1]]}) |> do.call(what = cbind)}

eigen_dec <- function(n, p, q) {

  D_mat <- build_D_mat(n, q)
  P  <- crossprod(D_mat)

  ei <- eigen(P, symmetric = TRUE)
  s <- sort(ei$values)[seq_len(p)]
  s[seq_len(q)] <- 0

  U <- ei$vectors[, order(ei$values)]

  X <- seq_len(q) |>
    lapply(\(k) (seq_len(n) - (n + 1) / 2) ^ (k - 1)) |>
    lapply(\(x) x / norm(x, "2")) |>
    do.call(what = cbind)
  Z <- U[, q + seq_len(p - q), drop = FALSE]

  list(X = X, Z = Z, s = s)
}

compute_res_deviance <- function(D, D_hat) {

  D_diff <- D - D_hat
  log_D_diff <- ifelse(D == 0, 0, D * (log(D) - log(D_hat)))

  sign(D_diff) * sqrt(2 * pmax(log_D_diff - D_diff, 0))
}

compute_deviance <- function(D, D_hat) {

  res_deviance <- compute_res_deviance(D, D_hat)

  sum(res_deviance * res_deviance)
}
