# Model matrices----

eigen_dec <- function(n, q) {

  P <- create_P_compact_cpp(n, q)
  s <- eigenvalues_compact_lapack(P)
  s[seq_len(q)] <- 0 # numerical issues without this line

  list(P = P, s = s)
}

get_model_matrices <- function(n, q) {

  if (length(n) == 1) {
    eig <- eigen_dec(n, q)
    P <- eig$P
    s <- eig$s
  } else {
    eig <- list(eigen_dec(n[[1]], q[[1]]), eigen_dec(n[[2]], q[[2]]))
    P <- list(eig[[1]]$P, eig[[2]]$P)
    s_eig <- list(eig[[1]]$s, eig[[2]]$s)
    s <- list(rep(s_eig[[1]], n[[2]]), rep(s_eig[[2]], each = n[[1]]))
  }
  list(P = P, s = s)
}

q_diff <- function(y_hat, q, margin = 0, final_transpose = TRUE) {

  if (!is.null(dim(y_hat)) && margin == 1) y_hat <- t(y_hat)
  out <- diff(y_hat, differences = q)
  if (final_transpose) out <- t(out)
  return(out)
}

get_RESS_2d_aux <- function(y_hat, q, margin) {

  delta <- q_diff(y_hat, q, margin, final_transpose = FALSE)
  sum(delta * delta)
}

get_RESS <- function(y_hat, q) {

  if (length(q) == 1) {
    delta <- q_diff(y_hat, q)
    sum(delta * delta)
  } else {
    c(get_RESS_2d_aux(y_hat, q[[1]], 2), get_RESS_2d_aux(y_hat, q[[2]], 1))
  }
}

get_edf_tot <- function(R, wt) {

  sum(wt * diag_V_compact_cpp(R))
}

# Model diagnoses----

compute_res_reg <- function(y, y_hat, wt) {

  y[wt == 0] <- 0
  sqrt(wt) * (y - y_hat)
}

#' Deviance residuals for Poisson GLM
#'
#' @param D Vector or matrix containing the number of observed events
#' @param D_hat Vector or matrix containing the number of predicted events
#'
#' @returns A vector or matrix (depending on the input type, will be a matrix
#'   if at least one of the input is) containing the deviance residuals
#' @keywords internal
#' @export
compute_res_deviance <- function(D, D_hat) {

  D_diff <- D - D_hat
  log_D_diff <- ifelse(D == 0, 0, D * (log(D) - log(D_hat)))

  sign(D_diff) * sqrt(2 * pmax(log_D_diff - D_diff, 0))
}

compute_K <- function(p, q, s, lambda, R) {

  if (all(p == q)) 0 else {

    log_det_S <- if (length(lambda) == 1) {

      (p - q) * log(lambda) + sum(log(s[- seq_len(q)]))

    } else {

      is_s_pos <- s[[1]] != 0 | s[[2]] != 0
      lambda_s <- list(lambda[[1]] * s[[1]], lambda[[2]] * s[[2]])
      (lambda_s[[1]] + lambda_s[[2]])[is_s_pos] |> log() |> sum()
    }
    diag_R <- R[nrow(R), ]
    2 * sum(log(diag_R)) - log_det_S - prod(q) * log(2 * pi)
  }
}

compute_REML <- function(dev, pen, p, q, s, lambda, R, l_sat = 0, K = NULL) {

  if (is.null(K)) K <- compute_K(p, q, s, lambda, R)
  (dev + pen) / 2 + K / 2 - l_sat
}

#' Diagnosis for Model Fit
#'
#' @param dev Deviance of the model
#' @param pen Penalization force
#' @param sum_edf Effective dimension
#' @param n_pos Number of strictly positive weights
#' @param p Number of parameters
#' @param q Orders of the penalization matrices
#' @param s eigenvalues for the penalization matrix
#' @param lambda smoothing parameters
#' @param R triangular factor
#' @param l_sat likelihood for the saturated model
#'
#' @returns A data.frame containing various diagnosis about the fit
#' @keywords internal
get_diagnosis <- function(dev, pen, sum_edf, n_pos, p, q, s, lambda, R, l_sat = 0) {

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2

  K <- compute_K(p, q, s, lambda, R)
  REML <- compute_REML(dev = dev, pen = pen, q = q, K = K, l_sat = l_sat)

  data.frame(dev = dev, pen = pen, K = K, sum_edf = sum_edf, n_pos = n_pos,
             AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)
}

# Update functions----

get_omega_pos <- function(s, lambda_s) {

  is_s_pos <- s[[1]] != 0 | s[[2]] != 0
  lambda_s_pos <- list(lambda_s[[1]][is_s_pos], lambda_s[[2]][is_s_pos])
  sum_lambda_s_pos <- lambda_s_pos[[1]] + lambda_s_pos[[2]]

  list(lambda_s_pos[[1]] / sum_lambda_s_pos, lambda_s_pos[[2]] / sum_lambda_s_pos)
}


# Misc functions----

backsolve2 <- function(R, y) {

  dim(y) <- NULL
  backsolve_compact_cpp(R, backsolve_compact_cpp(R, y, transpose = TRUE))
}

extend_vector <- function(v, dims) {

  if (length(v) < dims) c(v, rep(v[[length(v)]], dims - 1)) else v
}

init_y_hat <- function(d, off) {

  ifelse(off == - Inf, 0, log(pmax(d, 1e-8)) - off)
}

compute_z <- function(y_hat, new_wt, d) {

  y_hat + ifelse(new_wt == 0, 0, d / new_wt - 1)
}

get_l_sat <- function(n) {

  - n / 2 * log(2 * pi)
}

find_starting_lambda <- function(wt, diag_P, lambda_start_ratio) {

  if (is.list(diag_P)) {
    dim(wt) <- NULL
    c(find_starting_lambda(wt, diag_P[[1]], lambda_start_ratio),
      find_starting_lambda(wt, diag_P[[2]], lambda_start_ratio))
  } else {
    exp(stats::uniroot(\(x) {mean(wt / (wt + exp(x) * diag_P)) - lambda_start_ratio}, interval = c(- 6, 12) * log(10))$root)
  }
}

# Display----

print_fit_infos <- function(procedure, lambda, counter, metric_names, values, old_values = NULL,
                            valid_step = NULL, step_factor = NULL, final = FALSE, digits = 3) {

  if (final) {

    display_string <- paste0(
      procedure, " procedure completed in ", counter,
      " iterations, smoothing parameters: ",
      paste0(format(lambda, digits = digits), collapse = ", "))

  } else {

    display_string <- if (counter == 1) paste0("Starting ", procedure, " procedure\n") else ""
    display_string <- paste0(
      display_string, "Iteration ", counter, ", ",
      "smoothing parameters: ", paste0(format(lambda, digits = digits), collapse = ", "))
  }
  if (!is.null(valid_step)) display_string <- paste0(display_string, " (step-halving factor ", format(step_factor, digits = digits), ")")
  for (i in seq_along(metric_names)) {

    string_i <- if (counter == 1 || final) {

      paste0(if (final) "final" else "initial", " ", metric_names[[i]], ": ", format(values[[i]], digits = digits))

    } else {

      paste0(metric_names[[i]], ": ",
             format(old_values[[i]], digits = digits), " => ",
             format(values[[i]], digits = digits),
             if (is.null(valid_step)) "" else {if (valid_step) " (accepted)" else " (rejected)"})
    }
    display_string <- paste0(display_string, ", ", string_i)
  }
  display_string <- paste0(display_string, "\n")
  cat(display_string)
}

# Sparse matrix computation----

combine_lambda_P <- function(lambda, P) {

  if (length(lambda) == 1) return(lambda * P)
  lambda_P <- list(lambda[[1]] * P[[1]], lambda[[2]] * P[[2]])
  combine_P_compact(lambda_P)
}

combine_P_compact <- function(C) {

  p <- c(ncol(C[[1]]), ncol(C[[2]]))
  q <- c(nrow(C[[1]]) - 1, nrow(C[[2]]) - 1)

  q_total <- p[[1]] * q[[2]]
  C_out <- matrix(0, nrow = q_total + 1, ncol = p[[1]] * p[[2]])

  C_out[q_total + 1, ] <- c(C[[1]][q[[1]] + 1, 1:p[[1]]]) + rep(C[[2]][q[[2]] + 1,], each = p[[1]])

  for (k in seq_len(q[[1]])) {
    C_out[q_total - k + 1, ] <- c(rep(0, k), C[[1]][q[[1]] - k + 1, (k + 1):p[[1]]])
  }

  for (k in seq_len(q[[2]])) {
    C_out[q_total - k * p[[1]] + 1, ] <- rep(C[[2]][q[[2]] - k + 1,], each = p[[1]])
  }
  return(C_out)
}

get_y_new_compact <- function(y, lambda_P, R_22, n_tot, ind_fit) {

  y_aug <- numeric(prod(n_tot)); dim(y_aug) <- n_tot; y_aug[ind_fit] <- y
  y_aux <- get_prod_P_y_compact_cpp(y_aug, lambda_P)[- ind_fit]
  y_new <- - backsolve_compact_cpp(R_22, backsolve_compact_cpp(R_22, y_aux, transpose = TRUE))
  return(y_new)
}

get_diag_V_pred_compact <- function(R, lambda_P, R_22, n, n_tot, ind_fit) {

  prod_n <- prod(n)
  prod_n_tot <- prod(n_tot)
  diag_V_pred <- numeric(prod_n_tot - prod_n)
  for (j in 1:prod_n) {

    ej <- numeric(prod_n); ej[j] <- 1
    kj <- backsolve_compact_cpp(R, ej)
    aj <- get_y_new_compact(kj, lambda_P, R_22, n_tot, ind_fit)
    diag_V_pred <- diag_V_pred + aj * aj
  }
  return(diag_V_pred)
}

# Tests----

banded_to_compact <- function(P, q) {

  p <- nrow(P)
  C <- matrix(0, q + 1, p)

  for (k in 0:q) {
    C[q + 1 - k, (k + 1):p] <- P[seq.int(1 + p * k, p ^ 2 - k, by = p + 1)]
  }
  return(C)
}

compact_to_tri <- function(C) {

  p <- ncol(C)
  q <- nrow(C) - 1
  R <- matrix(0, p, p)

  for (k in 0:q) {
    R[seq.int(1 + p * k, p ^ 2 - k, by = p + 1)] <- C[q + 1 - k, (k + 1):p]
  }
  return(R)
}

compact_to_sym <- function(C) {

  p <- ncol(C)
  q <- nrow(C) - 1
  R <- matrix(0, p, p)

  for (k in 0:q) {
    R[seq.int(1 + p * k, p ^ 2 - k, by = p + 1)] <- R[seq.int(1 + k, p ^ 2 - k * p, by = p + 1)] <- C[q + 1 - k, (k + 1):p]
  }
  return(R)
}
