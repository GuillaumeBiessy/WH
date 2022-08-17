#' Compute Block Kronecket Product Matrix
#'
#' @param XZ A list whose components are list of 2 matrices
#'
#' @returns A matrix obtained by applying the Kronecker product to the 2 matrices
#'   contained in each element of \code{XZ} and then horizontally concatenating
#'   the result
compute_XZ_mat <- function(XZ) {
  lapply(XZ, \(X) {X[[2]] %x% X[[1]]}) |> do.call(what = cbind)}

#' Eigen Decomposition of Penalization Matrix
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
eigen_dec <- function(n, q, p) {

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

#' Deviance Residuals for Poisson GLM
#'
#' @param D Vector or matrix containing the number of observed events
#' @param D_hat Vector or matrix containing the number of predicted events
#'
#' @returns A vector or matrix (depending on the input type, will be a matrix
#'   if at least one of the input is) containing the deviance residuals
compute_res_deviance <- function(D, D_hat) {

  D_diff <- D - D_hat
  log_D_diff <- ifelse(D == 0, 0, D * (log(D) - log(D_hat)))

  sign(D_diff) * sqrt(2 * pmax(log_D_diff - D_diff, 0))
}

#' Deviance for Poisson GLM
#'
#' @inheritParams compute_res_deviance
#'
#' @returns The model deviance
compute_deviance <- function(D, D_hat) {

  res_deviance <- compute_res_deviance(D, D_hat)

  sum(res_deviance * res_deviance)
}
