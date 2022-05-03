build_D_q <- function(n, q) {diff(diag(n), differences = q)}

WH_1d_vintage <- function(y, wt = rep(1, length(y)), q = 2, lambda = 1e3, type_IC = "bayesian") {

  n <- length(y)
  W <- diag(wt) # weight matrix
  D_q <- build_D_q(n, q) # difference matrix
  P <- lambda * crossprod(D_q) # penalization matrix

  Psi <- solve(W + P) # variance / covariance matrix
  y_hat <- colSums(wt * y * Psi) # fitted value

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- sum(wt * diag(Psi)) # effective degrees of freedom
  sigma_2_hat <- sum(res ^ 2) / (n - edf) # estimated dispersion parameter
  std_y_hat <- switch(type_IC, # standard deviation of fit
                      freq = sqrt(sigma_2_hat * colSums(wt * t(Psi) * Psi)),
                      bayesian = sqrt(sigma_2_hat * diag(Psi)))

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(y) # set names for output vectors

  out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, lambda = lambda)
  class(out) <- "WH_1d"

  return(out)
}

WH_1d <- function(y, wt = rep(1, length(y)), q = 2, lambda_start = 1, n_iter = 30) {

  # Initialization
  n <- length(y)
  W <- diag(wt)
  D_q <- build_D_q(n, q) # difference matrix

  SVD <- SVD_aux(D_q, q)

  X <- SVD$X
  Z <- SVD$Z

  tXWX <- t(X) %*% (wt * X)
  tXWZ <- t(X) %*% (wt * Z)
  tZWZ <- t(Z) %*% (wt * Z)
  tXWy <- t(X) %*% (wt * y)
  tZWy <- t(Z) %*% (wt * y)

  dLambda <- diag(SVD$Sigma)

  lambda <- lambda_start

  # Loop
  for (i in seq(n_iter)) {

    Gm <- lambda * diag(dLambda)
    dG <- 1 / (dLambda * lambda)

    Cm <- solve(tZWZ + Gm)
    CmtB <- Cm %*% t(tXWZ)
    S <- solve(tXWX - tXWZ %*% CmtB)

    Psi_hat <- vector("list", 4) |>
      matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))

    Psi_hat[["X", "X"]] <- S
    Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
    Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
    Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]

    beta_hat  <- c(Psi_hat[["X", "X"]] %*% tXWy + Psi_hat[["X", "Z"]] %*% tZWy)
    alpha_hat <- c(Psi_hat[["Z", "X"]] %*% tXWy + Psi_hat[["Z", "Z"]] %*% tZWy)

    zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)

    ed_re <- lambda * sum(zeta * dG * dLambda)
    RRSS <- sum(dLambda * alpha_hat ^ 2)
    lambda <- ed_re / RRSS

    # cat("Lambda:", lambda, "\n")
  }

  y_hat <- c(X %*% beta_hat + Z %*% alpha_hat)

  Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
                   cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))

  XZ <- cbind(X, Z)
  std_y_hat <- sqrt(rowSums((XZ %*% Psi_hat) * XZ))

  edf <- q + ed_re
  res <- sqrt(wt) * (y - y_hat)

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(y) # set names for output vectors

  out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, lambda = lambda)
  class(out) <- "WH_1d"

  return(out)
}

plot.WH_1d <- function(x, trans) {

  if (missing(trans)) trans <- \(x) x

  plot(names(x$y), trans(x$y), xlab = "Age", ylab = "log - taux de décès")
  lines(names(x$y), trans(x$y_hat), col = "blue")
  lines(names(x$y), trans(x$y_hat - 2 * x$std_y_hat), col = "red", lty = 3)
  lines(names(x$y), trans(x$y_hat + 2 * x$std_y_hat), col = "red", lty = 3)
}

WH_2d_vintage <- function(y, wt, q = c(2, 2), lambda = c(1e2, 1e2)) {

  n <- dim(y)
  W <- diag(wt) # weight matrix
  D_q <- map2(n, q, build_D_q) # difference matrices
  P <- lambda[[1]] * diag(n[[2]]) %x% crossprod(D_q[[1]]) +
    lambda[[2]] * crossprod(D_q[[2]]) %x% diag(n[[1]]) # penalization matrix

  Psi <- solve(W + P) # variance / covariance matrix
  y_hat <- colSums(wt * y * Psi) # fitted value

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- sum(wt * diag(Psi)) # effective degrees of freedom
  sigma_2_hat <- sum(res ^ 2) / (n - edf) # estimated dispersion parameter
  std_y_hat <- switch(type_IC, # standard deviation of fit
                      freq = sqrt(sigma_2_hat * colSums(wt * t(Psi) * Psi)),
                      bayesian = sqrt(sigma_2_hat * diag(Psi)))

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(y) # set names for output vectors

  out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, lambda = lambda)
  class(out) <- "WH_2d"

  return(out)
}

WH_2d <- function(y, wt, q = c(2, 2),
                  lambda_start = c(1, 1), n_iter = 30) {

  # Initialization
  D_q <- map2(dim(y), q, build_D_q) # difference matrices
  SVD <- map2(D_q, q, SVD_aux)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  Sigma_SVD <- map(SVD, "Sigma")

  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]]))
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]]))

  G2XZ <- list(compute_G2XZ(X, X), compute_G2XZ(X, Z),
         compute_G2XZ(Z, X), compute_G2XZ(Z, Z)) |>
    matrix(2, 2, byrow = TRUE)
  tG2XZ <- G2XZ |> modify_depth(3, t)

  tXWX <- compute_tXWZ_GLAM(X, X, wt, tG2XZ[[1, 1]])
  tXWZ <- compute_tXWZ_GLAM(X, Z, wt, tG2XZ[[1, 2]])
  tZWZ <- compute_tXWZ_GLAM(Z, Z, wt, tG2XZ[[2, 2]])
  tXWz <- compute_tXWz_GLAM(X, y, wt)
  tZWz <- compute_tXWz_GLAM(Z, y, wt)

  var_comp <- build_var_comp_2d(X_SVD, Z_SVD, Sigma_SVD,
                                include_coef = list(Z = 3:5))

  dLambda <- var_comp$dLambda
  build_G <- var_comp$build_G

  lambda <- lambda_start

  # Loop
  for (i in seq(n_iter)) {

    G_built <- build_G(lambda)

    Gm <- G_built$Gm
    dG <- G_built$dG

    Cm <- Matrix::forceSymmetric(tZWZ + Gm) |>
      solve() |> as.matrix()
    CmtB <- Cm %*% t(tXWZ)
    S <- solve(tXWX - tXWZ %*% CmtB)

    Psi_hat <- vector("list", 4) |> matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))

    Psi_hat[["X", "X"]] <- S
    Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
    Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
    Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]

    beta_hat_vec  <- c(Psi_hat[["X", "X"]] %*% tXWz + Psi_hat[["X", "Z"]] %*% tZWz)
    alpha_hat_vec <- c(Psi_hat[["Z", "X"]] %*% tXWz + Psi_hat[["Z", "Z"]] %*% tZWz)

    zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)

    edr <- map2_dbl(dLambda, lambda, ~.y * sum(zeta * dG * .x))
    RRSS <- map_dbl(dLambda, ~sum(.x * alpha_hat_vec ^ 2))
    lambda <- edr / RRSS

    cat("Lambda:", lambda, "\n")
  }

  XZ <- c(X, Z)
  par_hat_vec <- c(beta_hat_vec, alpha_hat_vec)
  dims_par <- map(XZ, map_dbl, ncol)
  ind_par <- map_dbl(dims_par, prod) |> cumindex()

  beta_hat <- map2(ind_beta, dims_X, ~array(data = beta_hat_vec[.x], dim = .y))
  alpha_hat <- map2(ind_alpha, dims_Z, ~array(data = alpha_hat_vec[.x], dim = .y))

  par <- map2(ind_par, dims_par, ~array(data = par_hat_vec[.x], dim = .y))

  y_hat <- compute_Xtheta_GLAM(XZ, par)

  Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
                   cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))

  G2XZ <- compute_G2XZ(XZ, XZ)

  std_y_hat <- sqrt(compute_diagXPsitZ_GLAM(XZ, XZ, Psi_hat, G2XZ))

  y[D == 0] <- NA

  res <- sqrt(wt) * (y - y_hat)
  ed <- sum(q + edr)

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(y) # set names for output vectors
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(y) # set names for output vectors
  out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, ed = ed, lambda = lambda)
  class(out) <- "WH_2d"

  return(out)
}

plot.WH_2d <- function(x, trans) {

  if (missing(trans)) trans <- \(x) x

  contour(as.numeric(rownames(x$y_hat)),
          as.numeric(colnames(x$y_hat)),
          trans(x$y_hat))
}
