#' Build Difference Matrix of Given Order
#'
#' @inheritParams eigen_dec
#'
#' @returns The difference matrix of order \code{q} for a vector containing
#'   \code{n} observations.
#' @keywords internal
build_D_mat <- function(n, q) {diff(diag(n), differences = q)}

cum_index <- function(n) {

  map2(n, c(0, cumsum(n)[- length(n)]), \(x, y) y + seq_len(x))
}

#' Diagonal Stacking of Matrices
#'
#' @param ... One or several object of type `matrix`
#'
#' @returns A matrix obtained by stacking the input matrices diagonally
#' @keywords internal
blockdiag <- function(...) {

  L  <- Filter(Negate(is.null), list(...))

  n1 <- lapply(L, nrow)
  n2 <- lapply(L, ncol)

  M  <- matrix(0, do.call(sum, n1), do.call(sum, n2))

  i1 <- cum_index(n1)
  i2 <- cum_index(n2)

  for(i in seq_along(L)) M[i1[[i]], i2[[i]]] <- L[[i]]

  return(M)
}

#' Eigen Decomposition of Penalization Matrix
#'
#' @param n Number of observations in the problem
#' @param q Order of the penalization matrix
#' @param p Number of eigenvectors to keep
#'
#' @returns A list with components : - \code{U} a matrix whose columns are the
#'   eigenvectors of the penalization matrix - \code{s} a vector containing the
#'   eigenvalues sorted in ascending order
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

#' Deviance Residuals for Poisson GLM
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
#' @returns A data.frame containing various diagnoses about the fit
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

#' Lapply with Custom Return Type
#'
#' @param x A vector or list
#' @param f A function
#' @param output_type The desired return type. Should be one of "list" (the
#'   default), "integer", "numeric" or "character".
#'
#' @returns A list (if output_type = "list") or a vector, sharing names with x,
#'   built by applying f to each element of x. This is a rough implementation
#'   of the `map` function from the `purrr` package and aims at removing the
#'   dependency from that package.
#' @keywords internal
map <- function(x, f, output_type = "list") {

  nx <- length(x)
  l <- vector(output_type, nx)
  for (i in seq_len(nx)) {
    l[[i]] <- f(x[[i]])
  }
  names(l) <- names(x)

  return(l)
}

#' Bivariate Lapply with Custom Return Type
#'
#' @param x,y A couple of vectors or lists with the same size
#' @param f A function
#' @param output_type The desired return type. Should be one of "list" (the
#'   default), "integer", "numeric" or "character".
#'
#' @returns A list (if output_type = "list") or a vector, sharing names with x,
#'   built by applying f to each element of x and y simultaneously. This is a
#'   rough implementation of the `map2` function from the `purrr` package
#'   purrr aims at removing the dependency from that package.
#' @keywords internal
map2 <- function(x, y, f, output_type = "list") {

  nx <- length(x)
  ny <- length(y)

  if (nx != ny) stop("Inputs should have the same length")

  l <- vector(output_type, nx)
  for (i in seq_len(nx)) {
    l[[i]] <- f(x[[i]], y[[i]])
  }
  names(l) <- names(x)

  return(l)
}
