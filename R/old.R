# WH_1d_vintage <- function(D, E, q = 2, lambda = 1e3, method = "vintage", type_IC = "bayesian") {
#
#   y <- log(D / E) # observation vector
#   y[D == 0] <- 0
#   wt <- D
#
#   n <- length(y)
#   Dq <- diff(diag(n), differences = q)
#
#   switch(method,
#          vintage = {
#
#            S <- lambda * crossprod(Dq)
#            W <- diag(wt)
#            Psi <- solve(W + S)
#
#            y_hat <- colSums(wt * y * Psi)
#
#            edf <- sum(wt * diag(Psi))
#            res <- sqrt(wt) * (y - y_hat)
#            sigma_2_hat <- sum(res ^ 2) / (n - edf)
#            std_y_hat <- switch(type_IC,
#                                freq = sqrt(sigma_2_hat * colSums(wt * t(Psi) * Psi)),
#                                bayesian = sqrt(sigma_2_hat * diag(Psi)))
#          },
#          wood = {
#
#            sqwt <- sqrt(wt)
#            X_aug <- rbind(sqwt * diag(n), - sqrt(lambda) * Dq)
#
#            QR <- qr(X_aug)
#            Q <- qr.Q(QR)
#            R <- qr.R(QR)
#
#            f <- colSums((sqwt * y) * Q[seq_len(n),])
#
#            y_hat <- backsolve(R, f)
#            A <- tcrossprod(Q[seq_len(n),])
#
#            edf <- sum(diag(A))
#            res <- sqrt(wt) * (y - y_hat)
#            sigma_2_hat <- sum(res ^ 2) / (n - edf)
#            std_y_hat <- switch(type_IC,
#                                freq = sqrt(sigma_2_hat * colSums(t(A) * A) / wt),
#                                bayesian = sqrt(sigma_2_hat * diag(A) / wt))
#          })
#
#   names(y_hat) <- names(std_y_hat) <- names(res) <- names(y) # set names for output vectors
#
#   out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, lambda = lambda)
#   class(out) <- "WH_1d"
#
#   return(out)
# }
#
# WH_1d_performance <- function(D, E, q = 2, rho_start = 0, method = "vintage", type_IC = "bayesian") {
#
#   # Initialization
#   y <- log(D / E) # observation vector
#   y[D == 0] <- 0
#
#   n <- length(y)
#   wt <- D
#   sqwt <- sqrt(wt)
#
#   Dq <- diff(diag(n), differences = q)
#   S <- crossprod(Dq)
#
#   Q <- diag(n)
#   R <- sqwt * diag(n)
#
#   rho <- rho_start
#
#   # Loop
#
#   for(i in seq_len(2e2)) {
#
#     lambda <- exp(rho)
#     R_aug <- rbind(R, sqrt(lambda) * Dq)
#
#     SVD <- svd(R_aug)
#
#     U1 <- SVD$u[seq_len(n),]
#     Dm1 <- diag(1 / SVD$d)
#     V1 <- SVD$v
#
#     edf <- sum(U1 * U1)
#
#     y1 <- t(U1) %*% t(Q) %*% sqwt * y
#     M1 <- Dm1 %*% t(V1) %*% S %*% V1 %*% Dm1
#     K1 <- M1 %*% t(U1) %*% U1
#     K2 <- t(K1)
#
#     y_hat <- V1 %*% (diag(Dm1) * y1)
#     alpha <- sum(wt * (y - y_hat) ^ 2)
#     delta <- n - edf
#     sigma_2_hat <- alpha / delta
#     std_y_hat <- sqrt(sigma_2_hat * colSums(diag(Dm1) ^ 2 * t(V1) * t(V1)))
#
#     nu_g <- n * alpha / delta ^ 2
#     nu_u <- alpha / n - (2 * delta * sigma_2_hat) / n + sigma_2_hat
#     print(nu_g)
#
#     d_delta <- lambda * sum(diag(K1))
#     d2_delta <- d_delta - 2 * lambda ^ 2 * sum(t(M1) * K1)
#
#     d_alpha <- lambda * t(y1) %*% (2 * M1 - K1 - K2) %*% y1
#     d2_alpha <- 2 * lambda ^ 2 * t(y1) %*% (2 * M1 %*% K1 - 2 * M1 %*% M1 + 2 * K1 %*% M1) %*% y1
#
#     m_g <- n / delta ^ 2 * d_alpha - 2 * n * alpha / delta ^ 3 * d_delta
#     M_g <- - 2 * n / delta ^ 3 * d_delta %*% d_alpha + n / delta ^ 2 * d2_alpha - 2 * n / delta ^ 3 * d_alpha %*% d_delta +
#       6 * n * alpha / delta ^ 4 * d_delta %*% d_delta - 2 * n * alpha / delta ^ 3 * d2_delta
#
#     m_u <- d_alpha / n - 2 * sigma_2_hat / n * d_delta
#     M_u <- d2_alpha / n - 2 * sigma_2_hat / n * d2_delta
#
#     rho <- drop(rho - solve(M_g, m_g))
#   }
#
#   names(y_hat) <- names(std_y_hat) <- names(res) <- names(y) # set names for output vectors
#
#   out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, lambda = lambda)
#   class(out) <- "WH_1d"
#
#   return(out)
# }
#
# (bench::mark(WH_1d_vintage(D, E, method = "vintage"),
#              WH_1d_vintage(D, E, method = "wood")))
#
# WH_2d <- function(y, wt, q = c(2, 2), lambda_start = c(1, 1), n_iter = 10,
#                   type_IC = "bayesian", method = "GLAM") {
#
#   # Initialization
#   D_q <- map2(dim(y), q, build_D_q) # difference matrices
#   SVD <- map2(D_q, q, SVD_aux)
#
#   X_SVD <- map(SVD, "X")
#   Z_SVD <- map(SVD, "Z")
#   Sigma_SVD <- map(SVD, "Sigma")
#
#   X <- list(XX = list(X_SVD[[1]], X_SVD[[2]]))
#   Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
#             XZ = list(X_SVD[[1]], Z_SVD[[2]]),
#             ZZ = list(Z_SVD[[1]], Z_SVD[[2]]))
#
#   if (method == "GLM") {
#
#     X_mat <- compute_XZ_mat(X)
#     Z_mat <- compute_XZ_mat(Z)
#
#     tXWX <- compute_tXWZ_GLM(X_mat, X_mat, c(wt))
#     tXWZ <- compute_tXWZ_GLM(X_mat, Z_mat, c(wt))
#     tZWZ <- compute_tXWZ_GLM(Z_mat, Z_mat, c(wt))
#     tXWz <- compute_tXWz_GLM(X_mat, c(y), c(wt))
#     tZWz <- compute_tXWz_GLM(Z_mat, c(y), c(wt))
#
#     var_comp <- build_var_comp_2d(X_SVD, Z_SVD, Sigma_SVD,
#                                   include_coef = list(Z = 3:5))
#
#     dLambda <- var_comp$dLambda
#     build_G <- var_comp$build_G
#
#     lambda <- lambda_start
#
#     # Loop
#     for (i in seq(n_iter)) {
#
#       G_built <- build_G(lambda)
#
#       Gm <- G_built$Gm
#       dG <- G_built$dG
#
#       Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
#       CmtB <- Cm %*% t(tXWZ)
#       S <- solve(tXWX - tXWZ %*% CmtB)
#
#       Psi_hat <- vector("list", 4) |> matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))
#
#       Psi_hat[["X", "X"]] <- S
#       Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
#       Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
#       Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]
#
#       beta_hat_vec  <- c(Psi_hat[["X", "X"]] %*% tXWz + Psi_hat[["X", "Z"]] %*% tZWz)
#       alpha_hat_vec <- c(Psi_hat[["Z", "X"]] %*% tXWz + Psi_hat[["Z", "Z"]] %*% tZWz)
#
#       zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)
#
#       edr <- map2_dbl(dLambda, lambda, ~.y * sum(zeta * dG * .x))
#       RRSS <- map_dbl(dLambda, ~sum(.x * alpha_hat_vec ^ 2))
#       lambda <- edr / RRSS
#
#       cat("lambda:", round(lambda, 1), "\n")
#     }
#
#     XZ_mat <- cbind(X_mat, Z_mat)
#     par_hat_vec <- c(beta_hat_vec, alpha_hat_vec)
#
#     y_hat <- compute_Xtheta_GLM(XZ_mat, par_hat_vec)
#
#     Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
#                      cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))
#
#     std_y_hat <- sqrt(compute_diag_XPsitZ_GLM(XZ_mat, XZ_mat, Psi_hat))
#
#     res <- sqrt(wt) * (y - y_hat)
#     edf <- sum(q + edr)
#
#   } else {
#
#     G2XZ <- list(compute_G2XZ(X, X), compute_G2XZ(X, Z),
#                  compute_G2XZ(Z, X), compute_G2XZ(Z, Z)) |>
#       matrix(2, 2, byrow = TRUE)
#     tG2XZ <- G2XZ |> modify_depth(3, t)
#
#     tXWX <- compute_tXWZ_GLAM(X, X, wt, tG2XZ[[1, 1]])
#     tXWZ <- compute_tXWZ_GLAM(X, Z, wt, tG2XZ[[1, 2]])
#     tZWZ <- compute_tXWZ_GLAM(Z, Z, wt, tG2XZ[[2, 2]])
#     tXWz <- compute_tXWz_GLAM(X, y, wt)
#     tZWz <- compute_tXWz_GLAM(Z, y, wt)
#
#     var_comp <- build_var_comp_2d(X_SVD, Z_SVD, Sigma_SVD,
#                                   include_coef = list(Z = 3:5))
#
#     dLambda <- var_comp$dLambda
#     build_G <- var_comp$build_G
#
#     lambda <- lambda_start
#
#     QR <- qr(tZWZ)
#
#     # Loop
#     for (i in seq(n_iter)) {
#
#       G_built <- build_G(lambda)
#
#       Gm <- G_built$Gm
#       dG <- G_built$dG
#
#       Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
#       CmtB <- Cm %*% t(tXWZ)
#       S <- solve(tXWX - tXWZ %*% CmtB)
#
#       Psi_hat <- vector("list", 4) |> matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))
#
#       Psi_hat[["X", "X"]] <- S
#       Psi_hat[["X", "Z"]] <- - S %*% t(CmtB)
#       Psi_hat[["Z", "X"]] <- t(Psi_hat[["X", "Z"]])
#       Psi_hat[["Z", "Z"]] <- Cm - CmtB %*% Psi_hat[["X", "Z"]]
#
#       beta_hat_vec  <- c(Psi_hat[["X", "X"]] %*% tXWz + Psi_hat[["X", "Z"]] %*% tZWz)
#       alpha_hat_vec <- c(Psi_hat[["Z", "X"]] %*% tXWz + Psi_hat[["Z", "Z"]] %*% tZWz)
#
#       zeta <- colSums(t(Psi_hat[["Z", "X"]]) * tXWZ) + colSums(t(Psi_hat[["Z", "Z"]]) * tZWZ)
#
#       edr <- map2_dbl(dLambda, lambda, ~.y * sum(zeta * dG * .x))
#       RRSS <- map_dbl(dLambda, ~sum(.x * alpha_hat_vec ^ 2))
#       lambda <- edr / RRSS
#
#       cat("lambda:", round(lambda, 1), "\n")
#     }
#
#     XZ <- c(X, Z)
#
#     par_hat_vec <- c(beta_hat_vec, alpha_hat_vec)
#     dims_par <- map(XZ, map_dbl, ncol)
#     ind_par <- map_dbl(dims_par, prod) |> cumindex()
#     par <- map2(ind_par, dims_par, ~array(data = par_hat_vec[.x], dim = .y))
#
#     y_hat <- compute_Xtheta_GLAM(XZ, par)
#
#     Psi_hat <- rbind(cbind(Psi_hat[["X", "X"]], Psi_hat[["X", "Z"]]),
#                      cbind(Psi_hat[["Z", "X"]], Psi_hat[["Z", "Z"]]))
#
#     G2XZ <- compute_G2XZ(XZ, XZ)
#     std_y_hat <- sqrt(compute_diagXPsitZ_GLAM(XZ, XZ, Psi_hat, G2XZ))
#
#     res <- sqrt(wt) * (y - y_hat)
#     edf <- sum(q + edr)
#   }
#
#
#
#   dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(y)
#   dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(y) # set names for output vectors
#
#   out <- list(y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, lambda = lambda)
#   class(out) <- "WH_2d"
#
#   return(out)
# }
#
# # GLAM computation----
#
# vmix <- function(...) {c(t(do.call(list(...), what = cbind)))}
#
# G2 = function(X1, X2) {
#
#   (X1 %x% matrix(1, 1, ncol(X2))) * (matrix(1, 1, ncol(X1)) %x% X2)
# }
#
# Rotate = function(A) {
#
#   d <- seq_along(dim(A))
#   aperm(A, c(d[-1], d[1]))
# }
#
# tidyXM = function(XM, d) {
#
#   array(XM, c(nrow(XM), d[-1]))
# }
#
# Htransform = function(X, A) {
#
#   d = dim(A)
#   M = matrix(A, nrow = d[1])
#   if (nrow(X) == 1) {
#     if (all(X == 1)) return(tidyXM(t(colSums(M)), d))
#     else return(tidyXM(t(colSums(c(X) * M)), d))
#   }
#   if (nrow(X) == ncol(X)) {
#     X_min_diag <- X
#     diag(X_min_diag) <- 0
#     if (all(X_min_diag == 0)) {
#       Xdiag <- diag(X)
#       if (all(Xdiag == 1)) return(tidyXM(M, d))
#       else return(tidyXM(Xdiag * M, d))
#     }
#   }
#   return(tidyXM(X %*% M, d))
# }
#
# RH = function(X, A) {
#
#   Rotate(Htransform(X, A))
# }
#
# rec_RH = function(X, Theta) {
#
#   if (length(X) == 0) Theta else rec_RH(X[- 1], RH(X[[1]], Theta))
# }
#
# compute_G2XZ = function(X, Z) {
#
#   map(seq_along(X), function(i) {
#     map(seq_along(Z), function(j) {
#       map2(X[[i]], Z[[j]], G2)
#     })
#   }) |>
#     flatten() |>
#     matrix(nrow = length(X), byrow = TRUE)
# }
#
# compute_Xtheta_GLAM = function(X, Theta) {
#
#   c(map2(X, Theta, rec_RH) |> reduce(`+`))
# }
#
# compute_tXWz_GLAM = function(X, z, mu) {
#
#   X %>%
#     map_depth(2, t) %>%
#     map(rec_RH, z * mu) %>%
#     reduce(c) %>%
#     c()
# }
#
# compute_tXWZ_GLAM = function(X, Z, mu, tG2XZ) {
#
#   compute_tXWZ_GLAM_cell = function(Xi, Zj, mu, tG2_Xi_Zj) {
#
#     ll <- length(tG2_Xi_Zj)
#
#     pi <- map_dbl(Xi, ncol)
#     pj <- map_dbl(Zj, ncol)
#
#     candidate <- rec_RH(tG2_Xi_Zj, mu)
#     dim(candidate) <- vmix(pj, pi)
#     v <- c(seq(2, 2 * ll, 2), seq(1, 2 * ll - 1, 2))
#     candidate <- candidate %>% aperm(v)
#     dim(candidate) <- c(prod(pi), prod(pj))
#
#     return(candidate)
#   }
#
#   out <- map(seq_along(X), function(i) {
#     map(seq_along(Z), function(j) {
#       compute_tXWZ_GLAM_cell(X[[i]], Z[[j]], mu, tG2XZ[[i, j]])
#     }) %>% reduce(cbind)
#   }) %>% reduce(rbind)
#
#   return(out)
# }
#
# compute_diag_XPsitZ_GLAM = function(X, Z, PsiXZ, G2XZ) {
#
#   compute_diag_XPsitZ_GLAM_cell = function(Xi, Zj, Psi_Xi_Zj, G2_Xi_Zj) {
#
#     ll <- length(G2_Xi_Zj)
#
#     pi <- map_dbl(Xi, ncol)
#     pj <- map_dbl(Zj, ncol)
#
#     dim(Psi_Xi_Zj) <- c(pi, pj)
#     Psi_Xi_Zj <- aperm(Psi_Xi_Zj, vmix(ll + 1:ll, 1:ll))
#     dim(Psi_Xi_Zj) <- pi * pj
#
#     out <- rec_RH(G2_Xi_Zj, Psi_Xi_Zj)
#
#     return(out)
#   }
#
#   cumdX <- X %>% map(map_dbl, ncol) %>% map_dbl(prod) %>% cumindex
#   cumdZ <- Z %>% map(map_dbl, ncol) %>% map_dbl(prod) %>% cumindex
#
#   out <- map(seq_along(X), function(i) {
#     map(seq_along(Z), function(j) {
#       compute_diag_XPsitZ_GLAM_cell(X[[i]], Z[[j]], PsiXZ[cumdX[[i]], cumdZ[[j]]], G2XZ[[i, j]])
#     })
#   }) %>% flatten() %>% reduce(`+`) %>% c()
#
#   return(out)
# }
#
# A <- matrix(c(1,2,3,4,0,0), nrow = 2, byrow = TRUE)
# B <- matrix(c(5,6,7,0), nrow = 2, byrow = TRUE)
#
# `%r%` <- function(A, B) {
#
#   A[,rep(seq(ncol(A)), each = ncol(B))] * B[,rep(seq(ncol(B)), ncol(A))]
# }
#
# predict.WH_1d <- function(object, newdata = NULL, type_IC = "bayesian") {
#
#   data <- as.numeric(names(object$y))
#   full_data <- sort(union(data, newdata))
#   ind <- order(c(data, setdiff(full_data, data)))
#
#   n <- length(data)
#   n_pred <- length(full_data)
#   q <- object$q
#   lambda <- object$lambda
#
#   C1 <- diag(1, n, n_pred)[, ind] # constraint matrix
#   C2 <- matrix(0L, n_pred - n, n_pred)
#   C2[cbind(seq_len(n_pred - n), which(colSums(C1) == 0))] <- 1
#
#   if ("X" %in% names(object) && "Z" %in% names(object)) {
#
#     wt <- c(object$wt, rep(0, n_pred - n))
#
#     X <- object$X
#     Z <- object$Z
#
#     SVD_pred <- eigen_dec_pred(X, Z, n_pred, q, data, full_data)
#
#     D_mat_SVD <- SVD_pred$D_mat
#     X_pred <- SVD_pred$X
#     Z_pred <- SVD_pred$Z
#
#     Gm <- lambda * crossprod(D_mat_SVD %*% Z_pred)
#
#     tXWX <- compute_tXWZ(X_pred, X_pred, wt)
#     tXWZ <- compute_tXWZ(X_pred, Z_pred, wt)
#     tZWZ <- compute_tXWZ(Z_pred, Z_pred, wt)
#
#     Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
#     CmtB <- Cm %*% t(tXWZ)
#     S <- solve(tXWX - tXWZ %*% CmtB)
#
#     Psi_pred <- vector("list", 4) |>
#       matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))
#
#     Psi_pred[["X", "X"]] <- S
#     Psi_pred[["X", "Z"]] <- - S %*% t(CmtB)
#     Psi_pred[["Z", "X"]] <- t(Psi_pred[["X", "Z"]])
#     Psi_pred[["Z", "Z"]] <- Cm - CmtB %*% Psi_pred[["X", "Z"]]
#
#     Psi_pred <- rbind(cbind(Psi_pred[["X", "X"]], Psi_pred[["X", "Z"]]),
#                       cbind(Psi_pred[["Z", "X"]], Psi_pred[["Z", "Z"]]))
#
#     C1_XZ <- diag(1, ncol(X) + ncol(Z), ncol(X) + ncol(Z_pred))
#     C2_XZ <- matrix(0L, ncol(Z_pred) - ncol(Z), ncol(X) + ncol(Z_pred))
#     C2_XZ[cbind(seq_len(ncol(Z_pred) - ncol(Z)), which(colSums(C1_XZ) == 0))] <- 1
#
#     Psi_11_pred_inv <- (C1_XZ %*% Psi_pred %*% t(C1_XZ)) |> chol() |> chol2inv()
#     Psi_21_pred <- (C2_XZ %*% Psi_pred %*% t(C1_XZ))
#
#   } else {
#
#     wt <- t(C1) %*% c(object$wt)
#
#     W_pred <- diag(c(wt)) # extended weight matrix
#     D_mat_pred <- build_D_mat(n_pred, q) # extended difference matrices
#     P_pred <- lambda * crossprod(D_mat_pred) # extended penalization matrix
#
#     Psi_pred <- (W_pred + P_pred) |> chol() |> chol2inv() # un constrained variance / covariance matrix
#
#     Psi_11_pred_inv <- (C1 %*% Psi_pred %*% t(C1)) |> chol() |> chol2inv()
#     Psi_21_pred <- (C2 %*% Psi_pred %*% t(C1))
#   }
#
#   p <- length(object$par_hat)
#   XZ_pred <- rbind(diag(p), Psi_21_pred %*% Psi_11_pred_inv)
#
#   if ("X" %in% names(object) && "Z" %in% names(object)) XZ_pred <- cbind(X_pred, Z_pred) %*% XZ_pred
#
#   y_hat <- c(XZ_pred %*% object$par_hat)[ind]
#
#   sigma_2_hat <- if ("sigma_2_hat" %in% names(object)) object$sigma_2_hat else 1
#   std_y_hat <- sqrt(sigma_2_hat * compute_diag_XPsitZ(XZ_pred, XZ_pred, switch(
#     type_IC,
#     freq = object$Psi_hat %*% (c(object$wt) * object$Psi_hat),
#     bayesian = object$Psi_hat)))[ind]
#
#   names(y_hat) <- names(std_y_hat) <- full_data
#
#   object$y_hat <- y_hat
#   object$std_y_hat <- std_y_hat
#
#   return(object)
# }
#
# predict.WH_2d <- function(object, newdata = NULL, type_IC = "bayesian") {
#
#   data <- dimnames(object$y) |> map(as.numeric)
#   full_data <- map2(data, newdata, \(x,y) sort(union(x, y)))
#   ind <- map2(data, full_data, \(x,y) order(c(x, setdiff(y, x))))
#
#   n <- map_int(data, length)
#   n_pred <- map_int(full_data, length)
#   q <- object$q
#   lambda <- object$lambda
#
#   C1 <- map2(n, n_pred, \(x,y) diag(1, x, y)) |>
#     map2(ind, \(x,y) x[,y]) |>
#     rev() |>
#     reduce(kronecker) # constraint matrix
#   C2 <- matrix(0L, prod(n_pred) - prod(n), prod(n_pred))
#   C2[cbind(seq_len(prod(n_pred) - prod(n)), which(colSums(C1) == 0))] <- 1
#
#   if ("X" %in% names(object) && "Z" %in% names(object)) {
#
#     wt <- c(object$wt, rep(0, prod(n_pred) - prod(n)))
#
#     SVD_pred <- pmap(.l = list(X = object$X_SVD, Z = object$Z_SVD, n_pred = n_pred,
#                                q = q, data = data, full_data = full_data),
#                      .f = eigen_dec_pred)
#
#     D_mat_SVD <- map(SVD_pred, "D_mat")
#     X_SVD <- map(SVD_pred, "X")
#     Z_SVD <- map(SVD_pred, "Z")
#
#     X_pred <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
#       compute_XZ_mat()
#     Z_pred <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
#                    XZ = list(X_SVD[[1]], Z_SVD[[2]]),
#                    ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
#       compute_XZ_mat()
#
#     Gm <- build_var_comp_2d_pred(D_mat_SVD, X_SVD, Z_SVD, lambda)
#
#     tXWX <- compute_tXWZ(X_pred, X_pred, wt)
#     tXWZ <- compute_tXWZ(X_pred, Z_pred, wt)
#     tZWZ <- compute_tXWZ(Z_pred, Z_pred, wt)
#
#     Cm <- (tZWZ + Gm) |> chol() |> chol2inv()
#     CmtB <- Cm %*% t(tXWZ)
#     S <- solve(tXWX - tXWZ %*% CmtB)
#
#     Psi_pred <- vector("list", 4) |>
#       matrix(2, 2, dimnames = list("X", "Z") |> list() |> rep(2))
#
#     Psi_pred[["X", "X"]] <- S
#     Psi_pred[["X", "Z"]] <- - S %*% t(CmtB)
#     Psi_pred[["Z", "X"]] <- t(Psi_pred[["X", "Z"]])
#     Psi_pred[["Z", "Z"]] <- Cm - CmtB %*% Psi_pred[["X", "Z"]]
#
#     Psi_pred <- rbind(cbind(Psi_pred[["X", "X"]], Psi_pred[["X", "Z"]]),
#                       cbind(Psi_pred[["Z", "X"]], Psi_pred[["Z", "Z"]]))
#
#     build_C_block <- function(n, n_pred) {
#
#       map2(n, n_pred, \(x,y) diag(1, x, y)) %>% rev %>% reduce(`%x%`)
#     }
#
#     build_C <- function(n, n_pred) {
#
#       map2(n, n_pred, build_C_block) |> blockdiag()
#     }
#
#     XZ_dim <- list(XX = list(object$X_SVD[[1]], object$X_SVD[[2]]),
#                    ZX = list(object$Z_SVD[[1]], object$X_SVD[[2]]),
#                    XZ = list(object$X_SVD[[1]], object$Z_SVD[[2]]),
#                    ZZ = list(object$Z_SVD[[1]], object$Z_SVD[[2]])) |>
#       map(\(x) map_int(x, \(y) ncol(y)))
#
#     XZ_pred_dim <- list(XX = list(X_SVD[[1]], X_SVD[[2]]),
#                         ZX = list(Z_SVD[[1]], X_SVD[[2]]),
#                         XZ = list(X_SVD[[1]], Z_SVD[[2]]),
#                         ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
#       map(\(x) map_int(x, \(y) ncol(y)))
#
#     C1_XZ <- build_C(XZ_dim, XZ_pred_dim)
#     C2_XZ <- matrix(0L, ncol(Z_pred) - ncol(Z), ncol(C1_XZ))
#     C2_XZ[cbind(seq_len(ncol(Z_pred) - ncol(Z)), which(colSums(C1_XZ) == 0))] <- 1
#
#     Psi_11_pred_inv <- (C1_XZ %*% Psi_pred %*% t(C1_XZ)) |> chol() |> chol2inv()
#     Psi_21_pred <- (C2_XZ %*% Psi_pred %*% t(C1_XZ))
#
#     XZ_pred <- cbind(X_pred, Z_pred) %*% (t(C1_XZ) + t(C2_XZ) %*% Psi_21_pred %*% Psi_11_pred_inv)
#
#   } else {
#
#     wt <- t(C1) %*% c(object$wt)
#
#     W_pred <- diag(c(wt)) # extended weight matrix
#     D_mat_pred <- map2(n_pred, q, build_D_mat) # extended difference matrices
#     P_pred <- object$lambda[[1]] * diag(n_pred[[2]]) %x% crossprod(D_mat_pred[[1]]) +
#       object$lambda[[2]] * crossprod(D_mat_pred[[2]]) %x% diag(n_pred[[1]]) # extended penalization matrix
#
#     Psi_pred <- (W_pred + P_pred) |> chol() |> chol2inv() # un constrained variance / covariance matrix
#
#     Psi_11_pred_inv <- (C1 %*% Psi_pred %*% t(C1)) |> chol() |> chol2inv()
#     Psi_21_pred <- (C2 %*% Psi_pred %*% t(C1))
#
#     p <- length(object$par_hat)
#     XZ_pred <- rbind(diag(p), Psi_21_pred %*% Psi_11_pred_inv)
#   }
#
#   y_hat <- c(XZ_pred %*% object$par_hat)
#
#   sigma_2_hat <- if ("sigma_2_hat" %in% names(object)) object$sigma_2_hat else 1
#   std_y_hat <- sqrt(sigma_2_hat * compute_diag_XPsitZ(XZ_pred, XZ_pred, switch(
#     type_IC,
#     freq = object$Psi_hat %*% (c(object$wt) * object$Psi_hat),
#     bayesian = object$Psi_hat)))
#
#   y_hat_mat <- std_y_hat_mat <- rep(0, prod(n_pred))
#   dim(y_hat_mat) <- dim(std_y_hat_mat) <- map_int(full_data, length) # set dimension for output matrices
#   dimnames(y_hat_mat) <- dimnames(std_y_hat_mat) <- full_data # set names for output matrices
#
#   y_hat_mat[which(colSums(C1) == 1)] <- y_hat[seq_len(prod(n))]
#   y_hat_mat[which(colSums(C1) == 0)] <- y_hat[prod(n) + seq_len(prod(n_pred) - prod(n))]
#
#   std_y_hat_mat[which(colSums(C1) == 1)] <- std_y_hat[seq_len(prod(n))]
#   std_y_hat_mat[which(colSums(C1) == 0)] <- std_y_hat[prod(n) + seq_len(prod(n_pred) - prod(n))]
#
#   object$y_hat <- y_hat_mat
#   object$std_y_hat <- std_y_hat_mat
#
#   return(object)
# }
# eigen_dec_pred <- function(X, Z, n_pred, q, data, full_data) {
#
#   D_mat_pred <- build_D_mat(n_pred, q)
#
#   cols_fit <- which(full_data %in% data)
#   rows_fit <- cols_fit |> utils::head(- q)
#   cols_pred <- setdiff(seq_len(n_pred), cols_fit)
#   rows_pred <- setdiff(seq_len(n_pred - q), rows_fit)
#
#   D1 <- D_mat_pred[rows_pred, cols_fit, drop = FALSE]
#   D2 <- D_mat_pred[rows_pred, cols_pred, drop = FALSE]
#
#   D2_inv <- if (nrow(D2) > 0) solve(D2) else matrix(0, 0, 0)
#
#   D_mat_pred <- D_mat_pred[c(rows_fit, rows_pred), c(cols_fit, cols_pred)]
#
#   X_pred <- rbind(X, - D2_inv %*% D1 %*% X)
#   Z_pred <- blockdiag(list(Z, D2_inv))
#
#   list(D_mat = D_mat_pred, X = X_pred, Z = Z_pred)
# }
#
# build_var_comp_2d_pred <- function(D_mat, X, Z, lambda) {
#
#   F_blocks <- vector("list", 9)
#   dim(F_blocks) <- c(3, 3)
#
#   F_blocks <- F_blocks |> list() |> rep(2)
#
#   DZ  <- map2(D_mat, Z, `%*%`) |> map(crossprod)
#   tXX  <- map2(X, X, crossprod)
#   tZZ  <- map2(Z, Z, crossprod)
#   tXZ <- map2(X, Z, ~t(.x) %*% .y)
#
#   F_blocks[[1]][[1, 1]] <- tXX[[2]] %x% DZ[[1]]
#   F_blocks[[1]][[1, 3]] <- tXZ[[2]] %x% DZ[[1]]
#
#   F_blocks[[2]][[2, 2]] <- DZ[[2]] %x% tXX[[1]]
#   F_blocks[[2]][[2, 3]] <- DZ[[2]] %x% tXZ[[1]]
#
#   F_blocks[[1]][[3, 3]] <- tZZ[[2]] %x% DZ[[1]]
#   F_blocks[[2]][[3, 3]] <- DZ[[2]] %x% tZZ[[1]]
#
#   F_blocks[[1]][[3, 1]] <- t(F_blocks[[1]][[1, 3]])
#   F_blocks[[2]][[3, 2]] <- t(F_blocks[[2]][[2, 3]])
#
#   merge_matrix(F_blocks, lambda)
# }
#
# merge_matrix <- function(L, lambda) {
#
#   n_row <- L |> map(get_nrows) |> reduce(pmax)
#
#   L %>%
#     map(fill_L, n_row) %>%
#     map(crush_matrix) %>%
#     map2(lambda, `*`) %>%
#     reduce(`+`)
# }
#
# get_nrows <- function(L) {
#
#   L |> apply(1, compact) |> map_depth(2, nrow) |> map_dbl(reduce, max, .init = 0)
# }
#
# fill_L <- function(L, n_row) {
#
#   for(i in seq_along(n_row))
#     for(j in seq_along(n_row))
#       L[[i, j]] <- if (!is.null(L[[i, j]])) L[[i, j]] else matrix(0, n_row[[i]], n_row[[j]])
#
#   return(L)
# }
#
# crush_matrix <- function(ML) {
#
#   switch(ML |> dim() |> prod() |> as.character(),
#          "0" = NULL,
#          "1" = ML[[1, 1]],
#          ML |> apply(1, reduce, cbind) |> reduce(rbind))
# }
