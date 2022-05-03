cumindex <- function(n) {

  c(0, cumsum(n)[- length(n)]) |>
    map2(n, \(x,y) x + seq_len(y))
}

blockdiag <- function(L) {

  L  <- compact(L)
  n  <- length(L)

  n1 <- map(L, nrow)
  n2 <- map(L, ncol)

  i1 <- cumindex(n1)
  i2 <- cumindex(n2)

  M  <- matrix(0, do.call(sum, n1), do.call(sum, n2))
  for(i in seq_along(L)) M[i1[[i]], i2[[i]]] <- L[[i]]

  return(M)
}

SVD_aux <- function(D_q, q) {

  P  <- crossprod(D_q)
  p  <- nrow(P)
  ei <- eigen(P)
  U  <- ei$vectors[, order(ei$values)]

  X <- U[, seq_len(q), drop = FALSE]
  Z <- U[, q + seq_len(p - q), drop = FALSE]
  Sigma <- diag(sort(ei$values)[q + seq_len(p - q)])

  list(X = X, Z = Z, Sigma = Sigma)
}

# GLAM computation----

vmix <- function(...) {c(t(do.call(list(...), what = cbind)))}

G2 = function(X1, X2) {

  (X1 %x% matrix(1, 1, ncol(X2))) * (matrix(1, 1, ncol(X1)) %x% X2)
}

Rotate = function(A) {

  d <- seq_along(dim(A))
  aperm(A, c(d[-1], d[1]))
}

tidyXM = function(XM, d) {

  array(XM, c(nrow(XM), d[-1]))
}

Htransform = function(X, A) {

  d = dim(A)
  M = matrix(A, nrow = d[1])
  if (nrow(X) == 1) {
    if (all(X == 1)) return(tidyXM(t(colSums(M)), d))
    else return(tidyXM(t(colSums(c(X) * M)), d))
  }
  if (nrow(X) == ncol(X)) {
    X_min_diag <- X
    diag(X_min_diag) <- 0
    if (all(X_min_diag == 0)) {
      Xdiag <- diag(X)
      if (all(Xdiag == 1)) return(tidyXM(M, d))
      else return(tidyXM(Xdiag * M, d))
    }
  }
  return(tidyXM(X %*% M, d))
}

RH = function(X, A) {

  Rotate(Htransform(X, A))
}

rec_RH = function(X, Theta) {

  if (length(X) == 0) Theta else rec_RH(X[- 1], RH(X[[1]], Theta))
}

compute_G2XZ = function(X, Z) {

  map(seq_along(X), function(i) {
    map(seq_along(Z), function(j) {
      map2(X[[i]], Z[[j]], G2)
    })
  }) |>
    flatten() |>
    matrix(nrow = length(X), byrow = TRUE)
}

compute_Xtheta_GLAM = function(X, Theta) {

  c(map2(X, Theta, rec_RH) |> reduce(`+`))
}

compute_tXWz_GLAM = function(X, z, mu) {

  X %>%
    map_depth(2, t) %>%
    map(rec_RH, z * mu) %>%
    reduce(c) %>%
    c()
}

compute_tXWZ_GLAM = function(X, Z, mu, tG2XZ) {

  compute_tXWZ_GLAM_cell = function(Xi, Zj, mu, tG2_Xi_Zj) {

    ll <- length(tG2_Xi_Zj)

    pi <- map_dbl(Xi, ncol)
    pj <- map_dbl(Zj, ncol)

    candidate <- rec_RH(tG2_Xi_Zj, mu)
    dim(candidate) <- vmix(pj, pi)
    v <- c(seq(2, 2 * ll, 2), seq(1, 2 * ll - 1, 2))
    candidate <- candidate %>% aperm(v)
    dim(candidate) <- c(prod(pi), prod(pj))

    return(candidate)
  }

  out <- map(seq_along(X), function(i) {
    map(seq_along(Z), function(j) {
      compute_tXWZ_GLAM_cell(X[[i]], Z[[j]], mu, tG2XZ[[i, j]])
    }) %>% reduce(cbind)
  }) %>% reduce(rbind)

  return(out)
}

compute_diagXPsitZ_GLAM = function(X, Z, PsiXZ, G2XZ) {

  compute_diagXPsitZ_GLAM_cell = function(Xi, Zj, Psi_Xi_Zj, G2_Xi_Zj) {

    ll <- length(G2_Xi_Zj)

    pi <- map_dbl(Xi, ncol)
    pj <- map_dbl(Zj, ncol)

    dim(Psi_Xi_Zj) <- c(pi, pj)
    Psi_Xi_Zj <- aperm(Psi_Xi_Zj, vmix(ll + 1:ll, 1:ll))
    dim(Psi_Xi_Zj) <- pi * pj

    out <- rec_RH(G2_Xi_Zj, Psi_Xi_Zj)

    return(out)
  }

  cumdX <- X %>% map(map_dbl, ncol) %>% map_dbl(prod) %>% cumindex
  cumdZ <- Z %>% map(map_dbl, ncol) %>% map_dbl(prod) %>% cumindex

  out <- map(seq_along(X), function(i) {
    map(seq_along(Z), function(j) {
      compute_diagXPsitZ_GLAM_cell(X[[i]], Z[[j]], PsiXZ[cumdX[[i]], cumdZ[[j]]], G2XZ[[i, j]])
    })
  }) %>% flatten() %>% reduce(`+`) %>% c()

  return(out)
}

# Variance components----

build_var_comp_2d <- function(X, Z, Sigma, include_coef) {

  q <- map(X, ncol)
  p <- map(Z, ncol)

  # Marginal effects
  d1u <- Sigma[[1]]
  d2u <- Sigma[[2]]

  # Double interactions
  d12d <- diag(q[[2]]) %x% Sigma[[1]]
  d21d <- Sigma[[2]] %x% diag(q[[1]])

  d12dd <- diag(p[[2]]) %x% Sigma[[1]]
  d21dd <- Sigma[[2]] %x% diag(p[[1]])

  G_blocks <-
    list(d1u, d2u, d12d, d21d, d12dd) |>
    map(~0 * .x) |>
    list() |>
    rep(2)

  G_blocks[[1]][1]      <- list(d1u)
  G_blocks[[2]][2]      <- list(d2u)
  G_blocks[[1]][c(3,5)] <- list(d12d, d12dd)
  G_blocks[[2]][c(4,5)] <- list(d21d, d21dd)

  var_comp_builder(G_blocks, include_coef)
}

var_comp_builder <- function(G_blocks, include_coef) {

  G_blocks <- G_blocks |> map(~.x[include_coef$Z])

  dLambda <- G_blocks |>
    map_depth(2, diag) |>
    map(c, recursive = TRUE)

  build_G <- function(theta) {

    Gm <- G_blocks |>
      map2(theta, function(x, y) map(x, ~.x * y)) |>
      transpose() |>
      map(reduce, `+`) |>
      blockdiag()

    list(Gm = Gm, dG = 1 / diag(Gm))
  }

  list(dLambda = dLambda, build_G = build_G)
}
