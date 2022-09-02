#' Build difference matrix of order q
#'
#' @inheritParams eigen_dec
#'
#' @returns Difference matrix of order \code{q} for a vector containing \code{n}
#'   observations.
#' @keywords internal
build_D_mat <- function(n, q) {diff(diag(n), differences = q)}

#' Compute block Kronecker product matrix
#'
#' @param XZ A list whose components are list of 2 matrices
#'
#' @returns A matrix obtained by applying the Kronecker product to the 2 matrices
#'   contained in each element of \code{XZ} and then horizontally concatenating
#'   the result
#' @keywords internal
compute_XZ_mat <- function(XZ) {
  lapply(XZ, \(X) {X[[2]] %x% X[[1]]}) |> do.call(what = cbind)}

#' Eigen decomposition of penalization matrix
#'
#' @param n Number of observations in the problem
#' @param q Order of the penalization matrix
#' @param p Number of eigenvectors to keep
#'
#' @returns A list with components : - \code{X} a matrix whose columns are the
#'   eigenvectors associated with 0 eigenvalues - \code{Z} a matrix whose
#'   columns are the eigenvectors associated with non-0 eigenvalues sorted in
#'   ascending order - \code{s} a vector containing the eigenvalues sorted in
#'   ascending order
#' @keywords internal
eigen_dec <- function(n, q, p) {

  D_mat <- build_D_mat(n, q)
  P  <- crossprod(D_mat)

  ei <- eigen(P, symmetric = TRUE)
  s <- sort(ei$values)[seq_len(p)]
  s[seq_len(q)] <- 0

  X <- seq_len(q) |>
    lapply(\(k) (seq_len(n) - (n + 1) / 2) ^ (k - 1)) |>
    lapply(\(x) x / norm(x, "2")) |>
    do.call(what = cbind)
  Z <- ei$vectors[, order(ei$values)][, q + seq_len(p - q), drop = FALSE]

  out <- list(X = X, Z = Z, s = s)

  return(out)
}

cum_index <- function(n) {

  purrr::map2(n, c(0, cumsum(n)[- length(n)]), \(x, y) y + seq_len(x))
}

blockdiag <- function(...) {

  L  <- Filter(Negate(is.null), list(...))

  n1 <- purrr::map(L, nrow)
  n2 <- purrr::map(L, ncol)

  M  <- matrix(0, do.call(sum, n1), do.call(sum, n2))

  i1 <- cum_index(n1)
  i2 <- cum_index(n2)

  for(i in seq_along(L)) M[i1[[i]], i2[[i]]] <- L[[i]]

  return(M)
}

extend_eigen_dec <- function(data, full_data, q, p) {

  ind_fit <- which(full_data %in% data)
  ind_inf <- which(full_data < min(data))
  ind_sup <- which(full_data > max(data))

  n <- length(ind_fit)
  n_inf <- length(ind_inf)
  n_sup <- length(ind_sup)
  n_pred <- length(full_data)

  eig <- eigen_dec(n, q, p)

  D <- build_D_mat(n_pred, q)

  D1_inf <- D[ind_inf, ind_fit]
  D1_sup <- D[ind_sup - q, ind_fit]

  D2_inf <- D[ind_inf, ind_inf]
  D2_sup <- D[ind_sup - q, ind_sup]

  D2_inf_inv <- if (nrow(D2_inf) == 0) matrix(0, n_inf, n_inf) else solve(D2_inf)
  D2_sup_inv <- if (nrow(D2_sup) == 0) matrix(0, n_sup, n_sup) else solve(D2_sup)

  X <- matrix(0, n_pred, q)
  X[ind_fit, seq_len(q)] <- eig$X
  X[ind_inf, seq_len(q)] <- - D2_inf_inv %*% D1_inf %*% eig$X
  X[ind_sup, seq_len(q)] <- - D2_sup_inv %*% D1_sup %*% eig$X

  Z <- matrix(0, n_pred, p - q + n_inf + n_sup)
  Z[ind_fit, seq_len(p - q)] <- eig$Z
  Z[ind_inf, p - q + seq_len(n_inf)] <- D2_inf_inv
  Z[ind_sup, p - q + n_inf + seq_len(n_sup)] <- D2_sup_inv

  U <- cbind(X, Z)

  out <- list(U = U, D = D)

  return(out)
}

#' Deviance residuals for Poisson GLM
#'
#' @param D Vector or matrix containing the number of observed events
#' @param D_hat Vector or matrix containing the number of predicted events
#'
#' @returns A vector or matrix (depending on the input type, will be a matrix
#'   if at least one of the input is) containing the deviance residuals
#' @keywords internal
compute_res_deviance <- function(D, D_hat) {

  D_diff <- D - D_hat
  log_D_diff <- ifelse(D == 0, 0, D * (log(D) - log(D_hat)))

  out <- sign(D_diff) * sqrt(2 * pmax(log_D_diff - D_diff, 0))

  return(out)
}

#' Deviance for Poisson GLM
#'
#' @inheritParams compute_res_deviance
#'
#' @returns The model deviance
#' @keywords internal
compute_deviance <- function(D, D_hat) {

  res_deviance <- compute_res_deviance(D, D_hat)
  out <- sum(res_deviance * res_deviance)

  return(out)
}

#' Diagnosis for Model Fit
#'
#' @param dev Deviance of the model
#' @param pen Penalization force
#' @param sum_edf Effective dimension
#' @param n_pos Number of strictly positive weights
#' @param tr_log_P Trace of logarithm of penalization matrix
#' @param tr_log_Psi Trace of logarithm of variance-covariance matrix
#'
#' @returns A data.frame containing various diagnosis about the fit
#' @keywords internal
get_diagnosis <- function(dev, pen, sum_edf, n_pos, tr_log_P, tr_log_Psi) {

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- - (dev + pen - tr_log_P + tr_log_Psi) / 2

  out <- data.frame(dev = dev, pen = pen, sum_edf = sum_edf, n_pos = n_pos,
                    AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  return(out)
}

#' Diagnosis for Model Fit
#'
#' @param edf_par A vector containing effective degrees of freedom by parameter
#'   with fixed effects in front, fixed-random combinations, then random-random
#'   combinations
#' @param p Vector of number of parameters kept on each dimension
#' @param q Vector of penalization order on each dimension
#'
#' @returns A matrix containing effective degrees of freedom by parameter at the
#'   right position
#' @keywords internal
edf_par_to_matrix <- function(edf_par, p, q) {

  d <- p - q

  n_cum <- 0
  n_qq <- prod(q)
  mat_qq <- matrix(edf_par[seq_len(n_qq)], q[[1]], q[[2]])

  n_cum <- n_cum + n_qq
  n_dq <- d[[1]] * q[[2]]
  mat_dq <- matrix(edf_par[n_cum + seq_len(n_dq)], d[[1]], q[[2]])

  n_cum <- n_cum + n_dq
  n_qd <- q[[1]] * d[[2]]
  mat_qd <- matrix(edf_par[n_cum + seq_len(n_qd)], q[[1]], d[[2]])

  n_cum <- n_cum + n_qd
  n_dd <- d[[1]] * d[[2]]
  mat_dd <- matrix(edf_par[n_cum + seq_len(n_dd)], d[[1]], d[[2]])

  out <- cbind(rbind(mat_qq, mat_dq), rbind(mat_qd, mat_dd))

  return(out)
}
