#' Build difference matrix of order q
#'
#' @inheritParams eigen_dec
#'
#' @returns Difference matrix of order \code{q} for a vector containing \code{n}
#'   observations.
#' @keywords internal
build_D_mat <- function(n, q) {diff(diag(n), differences = q)}

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

  U <- ei$vectors[, order(ei$values)][, seq_len(p), drop = FALSE]
  U[, seq_len(q)] <- seq_len(q) |>
    lapply(\(k) (seq_len(n) - (n + 1) / 2) ^ (k - 1)) |>
    lapply(\(x) x / norm(x, "2")) |>
    do.call(what = cbind)

  s <- sort(ei$values)[seq_len(p)]
  s[seq_len(q)] <- 0

  out <- list(U = U, s = s)

  return(out)
}

extend_eigen_dec <- function(data, full_data, q, p, p_inf, p_sup) {

  ind_fit <- which(full_data %in% data)
  ind_inf <- which(full_data < min(data))
  ind_sup <- which(full_data > max(data))

  n <- length(ind_fit)
  n_inf <- length(ind_inf)
  n_sup <- length(ind_sup)
  n_pred <- length(full_data)

  if (missing(p_inf)) p_inf <- n_inf
  if (missing(p_sup)) p_sup <- n_sup

  eig <- eigen_dec(n, q, p)

  D <- build_D_mat(n_pred, q)

  D1_inf <- D[ind_inf, ind_fit]
  D1_sup <- D[ind_sup - q, ind_fit]

  D2_inf <- D[ind_inf, ind_inf]
  D2_sup <- D[ind_sup - q, ind_sup]

  D2_inf_inv <- if (nrow(D2_inf) == 0) matrix(0, n_inf, n_inf) else solve(D2_inf)
  D2_sup_inv <- if (nrow(D2_sup) == 0) matrix(0, n_sup, n_sup) else solve(D2_sup)

  U <- matrix(0, n_pred, p + p_inf + p_sup)
  U[ind_fit, seq_len(q)] <- eig$U[,seq_len(q)]
  U[ind_inf, seq_len(q)] <- - D2_inf_inv %*% D1_inf %*% eig$U[,seq_len(q)]
  U[ind_sup, seq_len(q)] <- - D2_sup_inv %*% D1_sup %*% eig$U[,seq_len(q)]
  U[ind_fit, q + seq_len(p - q)] <- eig$U[,q + seq_len(p - q)]
  U[ind_inf, p + seq_len(p_inf)] <- D2_inf_inv[, seq_len(p_inf)]
  U[ind_sup, p + p_inf + seq_len(p_sup)] <- D2_sup_inv[, seq_len(p_sup)]

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
