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

  out <- sign(D_diff) * sqrt(2 * round(log_D_diff - D_diff, 8))

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
