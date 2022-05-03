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
