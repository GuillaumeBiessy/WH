data(LTC)

library(Matrix)
library(purrr)

E <- LTC$E |> colSums()
D <- LTC$D |> colSums()

y <- D / E
y[D == 0] <- 0
n <- length(y)
w <- D
W <- diag(w)

lambda <- 1e2
porder <- 2
D <- diff(diag(n), differences = porder)
y_aug <- c(y, rep(0, n - porder))
X_aug <- rbind(diag(n), - lambda * D)
W_aug <- diag(c(w, rep(1, n - porder)))
QR <- qr(W_aug %*% X_aug)
Q <- qr.Q(QR)
R <- qr.R(QR)

f <- drop(t(Q) %*% (diag(W_aug) * y_aug))
y_hat <- drop(solve(R) %*% f)



WH_fit_1 <- WH_1d(D, E, lambda_start = 1)
WH_fit_2 <- WH_1d_vintage(D, E, lambda = WH_fit_1$lambda)
plot(WH_fit_1)
plot(WH_fit_2)

WH_fit_1 <- WH_2d(D, E)
WH_fit_2 <- WH_2d_vintage(D, E, lambda = WH_fit_1$lambda)

plot(WH_fit_1)
plot(WH_fit_2)

A <- matrix(c(1,2,3,4,0,0), nrow = 2, byrow = TRUE)
B <- matrix(c(5,6,7,0), nrow = 2, byrow = TRUE)

`%r%` <- function(A, B) {

  A[,rep(seq(ncol(A)), each = ncol(B))] * B[,rep(seq(ncol(B)), ncol(A))]
}
