# Regression----

WH_1d_reg_vintage <- function(y, wt = rep(1, length(y)),
                              lambda = 1e3, q = 2, type_IC = "bayesian") {

  n <- length(y)
  W <- diag(wt) # weight matrix
  D_mat <- build_D_mat(n, q) # difference matrix
  P <- lambda * crossprod(D_mat) # penalization matrix

  Psi_hat <- (W + P) |> chol() |> chol2inv()  # variance / covariance matrix
  par_hat <- y_hat <- c(Psi_hat %*% (wt * y)) # fitted value

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- sum(wt * diag(Psi_hat)) # effective degrees of freedom

  sigma_2_hat <- sum(res ^ 2) / (n - edf) # estimated dispersion parameter
  std_y_hat <- # standard deviation of fit
    sqrt(sigma_2_hat *
           switch(type_IC,
                  freq = rowSums(wt * Psi_hat * Psi_hat),
                  bayesian = diag(Psi_hat)))

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(wt) <- names(y) # set names for output vectors

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat, sigma_2_hat = sigma_2_hat,
              edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_1d"

  return(out)
}

WH_2d_reg_vintage <- function(y, wt = matrix(1, nrow = nrow(y), ncol = ncol(y)),
                              lambda = c(1e3, 1e3), q = c(2, 2), type_IC = "bayesian") {

  n <- dim(y)
  W <- diag(c(wt)) # weight matrix
  D_mat <- map2(n, q, build_D_mat) # difference matrices
  P <- lambda[[1]] * diag(n[[2]]) %x% crossprod(D_mat[[1]]) +
    lambda[[2]] * crossprod(D_mat[[2]]) %x% diag(n[[1]]) # penalization matrix

  Psi_hat <- (W + P) |> chol() |> chol2inv() # variance / covariance matrix
  par_hat <- y_hat <- c(Psi_hat %*% c(wt * y)) # fitted value

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- sum(wt * diag(Psi_hat)) # effective degrees of freedom

  sigma_2_hat <- sum(res ^ 2) / (prod(n) - edf) # estimated dispersion parameter
  std_y_hat <- # standard deviation of fit
    sqrt(sigma_2_hat *
           switch(type_IC,
                  freq = rowSums(c(wt) * Psi_hat * Psi_hat),
                  bayesian = diag(Psi_hat)))

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(wt) <- dim(y) # set dimension for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              Psi_hat = Psi_hat, sigma_2_hat = sigma_2_hat,
              par_hat = par_hat, edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_2d"

  return(out)
}

WH_1d_reg_mixed <- function(y, wt = rep(1, length(y)),
                            q = 2, type_IC = "bayesian") {

  # Initialization
  n <- length(y)
  SVD <- eigen_dec(n, q)

  X <- SVD$X
  Z <- SVD$Z
  dLambda <- diag(SVD$Sigma)

  XZ <- cbind(X, Z)

  tXWX <- compute_tXWZ(X, X, wt)
  tXWZ <- compute_tXWZ(X, Z, wt)
  tZWZ <- compute_tXWZ(Z, Z, wt)
  tXWy <- compute_tXWz(X, y, wt)
  tZWy <- compute_tXWz(Z, y, wt)

  lambda <- 1e2
  sigma_2_hat <- 1
  cond_lambda <- Inf

  # Loop
  while(cond_lambda > 1e-2) {

    Gm <- lambda * diag(dLambda)
    dG <- 1 / (dLambda * lambda)

    Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
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

    par_hat <- c(beta_hat, alpha_hat)
    y_hat <- c(XZ %*% par_hat)

    zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)
    ed_re <- lambda * sum(zeta * dG * dLambda)
    RRSS <- sum(dLambda * alpha_hat ^ 2)

    res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
    edf <- sum(q + ed_re) # effective degrees of freedom

    sigma_2_hat <- sigma_2_hat * sum(res ^ 2) / (n - edf)

    old_lambda <- lambda
    lambda <- ed_re / RRSS
    cond_lambda <- max(abs(log(lambda) - log(old_lambda)))
    print(paste("lambda:", format(lambda, digits = 2), "; sigma_2_hat:", format(sigma_2_hat, digits = 2)))
  }

  Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
                   cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))

  std_y_hat <- # standard deviation of fit
    sqrt(sigma_2_hat * switch(type_IC,
                              freq = compute_diag_XPsitZ(XZ, XZ, Psi_hat %*% (wt * Psi_hat)),
                              bayesian = compute_diag_XPsitZ(XZ, XZ, Psi_hat)))

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(wt) <- names(y) # set names for output vectors

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat, X = X, Z = Z,
              sigma_2_hat = sigma_2_hat, edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_1d"

  return(out)
}

WH_2d_reg_mixed <- function(y, wt = matrix(1, nrow = nrow(y), ncol = ncol(y)),
                            q = c(2, 2), type_IC = "bayesian") {

  # Initialization
  n <- dim(y)
  SVD <- map2(n, q, eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  Sigma_SVD <- map(SVD, "Sigma")

  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()

  tXWX <- compute_tXWZ(X, X, wt)
  tXWZ <- compute_tXWZ(X, Z, wt)
  tZWZ <- compute_tXWZ(Z, Z, wt)
  tXWy <- compute_tXWz(X, c(y), c(wt))
  tZWy <- compute_tXWz(Z, c(y), c(wt))

  var_comp <- build_var_comp_2d(X_SVD, Z_SVD, Sigma_SVD)
  dLambda <- var_comp$dLambda
  build_G <- var_comp$build_G

  lambda <- c(1e3, 1e3)
  sigma_2_hat <- 1
  cond_lambda <- Inf

  # Loop
  while(cond_lambda > 1e-2) {

    G_built <- build_G(lambda)

    Gm <- G_built$Gm
    dG <- G_built$dG

    Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
    CmtB <- Cm %*% t(tXWZ)
    S <- solve(tXWX - tXWZ %*% CmtB)

    Psi_hat <- vector("list", 4) |> matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))

    Psi_hat[["X", "X"]] <- S
    Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
    Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
    Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]

    beta_hat  <- c(Psi_hat[["X", "X"]] %*% tXWy + Psi_hat[["X", "Z"]] %*% tZWy)
    alpha_hat <- c(Psi_hat[["Z", "X"]] %*% tXWy + Psi_hat[["Z", "Z"]] %*% tZWy)

    y_hat <- c(X %*% beta_hat + Z %*% alpha_hat)

    zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)
    ed_re <- map2_dbl(dLambda, lambda, ~.y * sum(zeta * dG * .x))
    RRSS <- map_dbl(dLambda, ~sum(.x * alpha_hat * alpha_hat))

    res <- sqrt(c(wt)) * (y - y_hat) # (weighted) residuals
    edf <- sum(q + ed_re) # effective degrees of freedom

    sigma_2_hat <- sigma_2_hat * sum(res ^ 2) / (prod(n) - edf)

    old_lambda <- lambda
    lambda <- ed_re / RRSS
    cond_lambda <- max(abs(log(lambda) - log(old_lambda)))
    print(paste0(paste("lambda:", format(lambda, digits = 2)),
                 paste("; sigma_2_hat:", format(sigma_2_hat, digits = 2))))
  }

  XZ <- cbind(X, Z)
  par_hat <- c(beta_hat, alpha_hat)
  Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
                   cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))
  std_y_hat <- # standard deviation of fit
    sqrt(sigma_2_hat * switch(type_IC,
                              freq = compute_diag_XPsitZ(XZ, XZ, Psi_hat %*% (c(wt) * Psi_hat)),
                              bayesian = compute_diag_XPsitZ(XZ, XZ, Psi_hat)))

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(wt) <- dim(y) # set dimension for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat, X = X, Z = Z,
              sigma_2_hat = sigma_2_hat, edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_2d"

  return(out)
}

# Survival analysis----

WH_1d_surv_vintage <- function(ec, d, lambda = 1e3, q = 2, type_IC = "bayesian") {

  # Initialization
  n <- length(d)
  off <- log(ec |> pmax(1e-4))
  y <- ifelse(d == 0, NA, log(d) - off)
  z <- log(d |> pmax(1e-8)) - off
  wt <- d
  W <- diag(wt)

  D_mat <- build_D_mat(n, q) # difference matrix
  P <- lambda * crossprod(D_mat) # penalization matrix

  deviance <- cond_deviance <- Inf

  # Loop
  while(cond_deviance > 1e-4) {

    Psi_hat <- (W + P) |> chol() |> chol2inv()  # variance / covariance matrix
    par_hat <- y_hat <- c(Psi_hat %*% (wt * z)) # fitted value
    wt <- exp(y_hat + off)
    W <- diag(wt)
    z <- y_hat + d / wt - 1

    old_deviance <- deviance
    deviance <- compute_deviance(d, wt)
    cond_deviance <- abs(old_deviance / deviance - 1)
    print(paste("deviance:", format(deviance, digits = 3)))
  }

  res <- compute_res_deviance(d, wt) # (weighted) residuals
  edf <- sum(wt * diag(Psi_hat)) # effective degrees of freedom

  std_y_hat <- # standard deviation of fit
    sqrt(switch(type_IC,
                freq = rowSums(wt * Psi_hat * Psi_hat),
                bayesian = diag(Psi_hat)))

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(y) # set names for output vectors

  out <- list(y = y, z = z, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat,
              edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_1d"

  return(out)
}

WH_2d_surv_vintage <- function(ec, d, lambda = c(1e3, 1e3), q = c(2, 2), type_IC = "bayesian") {

  # Initialization
  n <- dim(d)
  off <- log(ec |> pmax(1e-4))
  y <- ifelse(d == 0, NA, log(d) - off)
  z <- log(d |> pmax(1e-8)) - off
  wt <- d

  W <- diag(c(wt))

  D_mat <- map2(n, q, build_D_mat) # difference matrices
  P <- lambda[[1]] * diag(n[[2]]) %x% crossprod(D_mat[[1]]) +
    lambda[[2]] * crossprod(D_mat[[2]]) %x% diag(n[[1]]) # penalization matrix

  deviance <- cond_deviance <- Inf

  # Loop
  while(cond_deviance > 1e-4) {

    Psi_hat <- (W + P) |> chol() |> chol2inv()  # variance / covariance matrix
    par_hat <- y_hat <- c(Psi_hat %*% c(wt * z)) # fitted value
    wt <- exp(y_hat + off)
    W <- diag(c(wt))
    z <- y_hat + d / wt - 1

    old_deviance <- deviance
    deviance <- compute_deviance(d, wt)
    cond_deviance <- abs(old_deviance / deviance - 1)
    print(paste("deviance:", format(deviance, digits = 2)))
  }

  res <- compute_res_deviance(d, wt) # (weighted) residuals
  edf <- sum(c(wt) * diag(Psi_hat)) # effective degrees of freedom

  std_y_hat <- # standard deviation of fit
    sqrt(switch(type_IC,
                freq = rowSums(c(wt) * Psi_hat * Psi_hat),
                bayesian = diag(Psi_hat)))

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <-
    dim(wt) <- dim(z) <- dim(y) # set dimension for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <-
    dimnames(wt) <- dimnames(z) <- dimnames(y) # set names for output matrices

  out <- list(y = y, z = z, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat,
              edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_2d"

  return(out)
}

WH_1d_surv_mixed <- function(ec, d, q = 2, type_IC = "bayesian") {

  # Initialization
  n <- length(d)
  off <- log(ec |> pmax(1e-4))
  y <- ifelse(d == 0, NA, log(d) - off)
  z <- log(d |> pmax(1e-8)) - off
  wt <- d

  SVD <- eigen_dec(n, q)

  X <- SVD$X
  Z <- SVD$Z
  dLambda <- diag(SVD$Sigma)

  lambda <- 1e2
  deviance <- cond_deviance <- Inf

  # Loop
  while(cond_deviance > 1e-4) {

    tXWX <- compute_tXWZ(X, X, wt)
    tXWZ <- compute_tXWZ(X, Z, wt)
    tZWZ <- compute_tXWZ(Z, Z, wt)
    tXWz <- compute_tXWz(X, z, wt)
    tZWz <- compute_tXWz(Z, z, wt)

    cond_lambda <- Inf

    # Loop
    while(cond_lambda > 1e-2) {

      Gm <- lambda * diag(dLambda)
      dG <- 1 / (dLambda * lambda)

      Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
      CmtB <- Cm %*% t(tXWZ)
      S <- solve(tXWX - tXWZ %*% CmtB)

      Psi_hat <- vector("list", 4) |>
        matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))

      Psi_hat[["X", "X"]] <- S
      Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
      Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
      Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]

      beta_hat  <- c(Psi_hat[["X", "X"]] %*% tXWz + Psi_hat[["X", "Z"]] %*% tZWz)
      alpha_hat <- c(Psi_hat[["Z", "X"]] %*% tXWz + Psi_hat[["Z", "Z"]] %*% tZWz)

      zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)

      ed_re <- lambda * sum(zeta * dG * dLambda)
      RRSS <- sum(dLambda * alpha_hat * alpha_hat)

      old_lambda <- lambda
      lambda <- ed_re / RRSS
      cond_lambda <- max(abs(log(lambda) - log(old_lambda)))
      print(paste("lambda:", format(lambda, digits = 2)))
    }

    y_hat <- c(X %*% beta_hat + Z %*% alpha_hat) # fitted value
    wt <- exp(y_hat + off)
    z <- y_hat + d / wt - 1

    old_deviance <- deviance
    deviance <- compute_deviance(d, wt)
    cond_deviance <- abs(old_deviance / deviance - 1)
    print(paste("deviance:", format(deviance, digits = 2)))
  }

  XZ <- cbind(X, Z)
  par_hat <- c(beta_hat, alpha_hat)
  Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
                   cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))

  std_y_hat <- # standard deviation of fit
    sqrt(switch(type_IC,
                freq = compute_diag_XPsitZ(XZ, XZ, Psi_hat %*% (wt * Psi_hat)),
                bayesian = compute_diag_XPsitZ(XZ, XZ, Psi_hat)))

  res <- compute_res_deviance(d, wt) # (weighted) residuals
  edf <- sum(q + ed_re) # effective degrees of freedom

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(wt) <- names(y) # set names for output vectors

  out <- list(y = y, z = z, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat, X = X, Z = Z,
              edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_1d"

  return(out)
}

WH_2d_surv_mixed <- function(ec, d, q = c(2, 2), type_IC = "bayesian") {

  # Initialization
  n <- dim(d)
  off <- log(ec |> pmax(1e-4))
  y <- ifelse(d == 0, NA, log(d) - off)
  z <- log(d |> pmax(1e-8)) - off
  wt <- d

  SVD <- map2(n, q, eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  Sigma_SVD <- map(SVD, "Sigma")

  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()

  var_comp <- build_var_comp_2d(X_SVD, Z_SVD, Sigma_SVD)
  dLambda <- var_comp$dLambda
  build_G <- var_comp$build_G

  lambda <- c(1e3, 1e3)
  deviance <- cond_deviance <- Inf

  # Loop
  while(cond_deviance > 1e-4) {

    tXWX <- compute_tXWZ(X, X, wt)
    tXWZ <- compute_tXWZ(X, Z, wt)
    tZWZ <- compute_tXWZ(Z, Z, wt)
    tXWz <- compute_tXWz(X, z, wt)
    tZWz <- compute_tXWz(Z, z, wt)

    cond_lambda <- Inf

    # Loop
    while(cond_lambda > 1e-2) {

      G_built <- build_G(lambda)
      Gm <- G_built$Gm
      dG <- G_built$dG

      Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
      CmtB <- Cm %*% t(tXWZ)
      S <- solve(tXWX - tXWZ %*% CmtB)

      Psi_hat <- vector("list", 4) |>
        matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))

      Psi_hat[["X", "X"]] <- S
      Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
      Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
      Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]

      beta_hat  <- c(Psi_hat[["X", "X"]] %*% tXWz + Psi_hat[["X", "Z"]] %*% tZWz)
      alpha_hat <- c(Psi_hat[["Z", "X"]] %*% tXWz + Psi_hat[["Z", "Z"]] %*% tZWz)

      zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)

      ed_re <- map2_dbl(dLambda, lambda, ~.y * sum(zeta * dG * .x))
      RRSS <- map_dbl(dLambda, ~sum(.x * alpha_hat * alpha_hat))

      old_lambda <- lambda
      lambda <- ed_re / RRSS
      cond_lambda <- max(abs(log(lambda) - log(old_lambda)))
      print(paste("lambda:", format(lambda, digits = 2)))
    }

    y_hat <- c(X %*% beta_hat + Z %*% alpha_hat) # fitted value
    wt <- exp(y_hat + off)
    z <- y_hat + d / wt - 1

    old_deviance <- deviance
    deviance <- compute_deviance(d, wt)
    cond_deviance <- abs(old_deviance / deviance - 1)
    print(paste("deviance:", format(deviance, digits = 2)))
  }

  XZ <- cbind(X, Z)
  par_hat <- c(beta_hat, alpha_hat)
  Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
                   cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))

  std_y_hat <- # standard deviation of fit
    sqrt(switch(type_IC,
                freq = compute_diag_XPsitZ(XZ, XZ, Psi_hat %*% (c(wt) * Psi_hat)),
                bayesian = compute_diag_XPsitZ(XZ, XZ, Psi_hat)))

  res <- compute_res_deviance(d, wt) # (weighted) residuals
  edf <- sum(q + ed_re) # effective degrees of freedom

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <-
    dim(wt) <- dim(z) <- dim(y) # set dimension for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <-
    dimnames(wt) <- dimnames(z) <- dimnames(y) # set names for output matrices

  out <- list(y = y, z = z, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res,
              par_hat = par_hat, Psi_hat = Psi_hat, X = X, Z = Z,
              edf = edf, lambda = lambda, q = q)
  class(out) <- "WH_2d"

  return(out)
}

# Plots----

plot.WH_1d <- function(x, trans) {

  if (missing(trans)) trans <- \(x) x

  plot(names(x$y), trans(x$y), xlab = "Age", ylab = "log - taux de décès",
       xlim = range(as.numeric(names(x$y_hat))),
       ylim = trans(range(c(x$y_hat - 2 * x$std_y_hat),
                          c(x$y_hat + 2 * x$std_y_hat))))
  lines(names(x$y_hat), trans(x$y_hat), col = "blue")
  lines(names(x$y_hat), trans(x$y_hat - 2 * x$std_y_hat), col = "red", lty = 3)
  lines(names(x$y_hat), trans(x$y_hat + 2 * x$std_y_hat), col = "red", lty = 3)
}

plot.WH_2d <- function(x, trans) {

  if (missing(trans)) trans <- \(x) x

  contour(as.numeric(colnames(x$y_hat)),
          as.numeric(rownames(x$y_hat)),
          nlevels = 20,
          trans(t(x$y_hat)))
}

# Extrapolation----

predict.WH_1d <- function(object, newdata = NULL, type_IC = "bayesian") {

  data <- as.numeric(names(object$y))
  full_data <- sort(union(data, newdata))
  ind <- order(c(data, setdiff(full_data, data)))

  n <- length(data)
  n_pred <- length(full_data)

  C1 <- diag(1, n, n_pred)[, ind] # constraint matrix
  C2 <- matrix(0L, n_pred - n, n_pred)
  C2[cbind(seq_len(n_pred - n), which(colSums(C1) == 0))] <- 1

  wt <- t(C1) %*% c(object$wt)
  W_pred <- diag(c(wt)) # extended weight matrix

  D_mat_pred <- build_D_mat(n_pred, object$q) # extended difference matrices
  P_pred <- object$lambda * crossprod(D_mat_pred) # extended penalization matrix

  Psi_pred <- (W_pred + P_pred) |> chol() |> chol2inv() # un constrained variance / covariance matrix

  Psi_11_pred_inv <- (C1 %*% Psi_pred %*% t(C1)) |> chol() |> chol2inv()
  Psi_21_pred <- (C2 %*% Psi_pred %*% t(C1))

  p <- length(object$par_hat)
  XZ_pred <- rbind(diag(p), Psi_21_pred %*% Psi_11_pred_inv)
  if ("X" %in% names(object) && "Z" %in% names(object)) XZ_pred <- XZ_pred %*% cbind(object$X, object$Z)

  y_hat <- c(XZ_pred %*% object$par_hat)[ind]

  sigma_2_hat <- if ("sigma_2_hat" %in% names(object)) object$sigma_2_hat else 1
  std_y_hat <- sqrt(sigma_2_hat * compute_diag_XPsitZ(XZ_pred, XZ_pred, switch(
    type_IC,
    freq = object$Psi_hat %*% (c(object$wt) * object$Psi_hat),
    bayesian = object$Psi_hat)))[ind]

  names(y_hat) <- names(std_y_hat) <- full_data

  object$y_hat <- y_hat
  object$std_y_hat <- std_y_hat

  return(object)
}

predict.WH_2d <- function(object, newdata = NULL, type_IC = "bayesian") {

  data <- dimnames(object$y) |> map(as.numeric)
  full_data <- map2(data, newdata, \(x,y) sort(union(x, y)))
  ind <- map2(data, full_data, \(x,y) order(c(x, setdiff(y, x))))

  n <- map_int(data, length)
  n_pred <- map_int(full_data, length)

  C1 <- map2(n, n_pred, \(x,y) diag(1, x, y)) |>
    map2(ind, \(x,y) x[,y]) |>
    rev() |>
    reduce(kronecker) # constraint matrix
  C2 <- matrix(0L, prod(n_pred) - prod(n), prod(n_pred))
  C2[cbind(seq_len(prod(n_pred) - prod(n)), which(colSums(C1) == 0))] <- 1

  wt <- t(C1) %*% c(object$wt)
  W_pred <- diag(c(wt)) # extended weight matrix

  D_mat_pred <- map2(n_pred, object$q, build_D_mat) # extended difference matrices
  P_pred <- object$lambda[[1]] * diag(n_pred[[2]]) %x% crossprod(D_mat_pred[[1]]) +
    object$lambda[[2]] * crossprod(D_mat_pred[[2]]) %x% diag(n_pred[[1]]) # extended penalization matrix

  Psi_pred <- (W_pred + P_pred) |> chol() |> chol2inv() # un constrained variance / covariance matrix

  Psi_11_pred_inv <- (C1 %*% Psi_pred %*% t(C1)) |> chol() |> chol2inv()
  Psi_21_pred <- (C2 %*% Psi_pred %*% t(C1))

  p <- length(object$par_hat)

  XZ_pred <- rbind(diag(p), Psi_21_pred %*% Psi_11_pred_inv)
  if ("X" %in% names(object) && "Z" %in% names(object)) XZ_pred <- XZ_pred %*% cbind(object$X, object$Z)

  y_hat <- c(XZ_pred %*% object$par_hat)

  sigma_2_hat <- if ("sigma_2_hat" %in% names(object)) object$sigma_2_hat else 1
  std_y_hat <- sqrt(sigma_2_hat * compute_diag_XPsitZ(XZ_pred, XZ_pred, switch(
    type_IC,
    freq = object$Psi_hat %*% (c(object$wt) * object$Psi_hat),
    bayesian = object$Psi_hat)))

  y_hat_mat <- std_y_hat_mat <- rep(0, prod(n_pred))
  dim(y_hat_mat) <- dim(std_y_hat_mat) <- map_int(full_data, length) # set dimension for output matrices
  dimnames(y_hat_mat) <- dimnames(std_y_hat_mat) <- full_data # set names for output matrices

  y_hat_mat[which(colSums(C1) == 1)] <- y_hat[seq_len(prod(n))]
  y_hat_mat[which(colSums(C1) == 0)] <- y_hat[prod(n) + seq_len(prod(n_pred) - prod(n))]

  std_y_hat_mat[which(colSums(C1) == 1)] <- std_y_hat[seq_len(prod(n))]
  std_y_hat_mat[which(colSums(C1) == 0)] <- std_y_hat[prod(n) + seq_len(prod(n_pred) - prod(n))]

  object$y_hat <- y_hat_mat
  object$std_y_hat <- std_y_hat_mat

  return(object)
}
