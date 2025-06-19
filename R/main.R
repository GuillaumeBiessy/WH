# Wrappers functions----

#' Whittaker-Henderson Smoothing
#'
#' Main package function to apply Whittaker-Henderson Smoothing in a survival
#' analysis framework. It takes as input two vectors / matrices of observed
#' events and associated central exposure and estimate a smooth version of the
#' log-hazard rate. Smoothing parameters may be supplied or automatically chosen
#' according to a specific criterion such as `"REML"` (recommended), `"AIC"`,
#' `"BIC"` or `"GCV"`. Whittaker-Henderson Smoothing may be applied in a full
#' maximum likelihood framework (strongly recommended) or an asymptotic
#' (approximate) Gaussian framework.
#'
#' @param d Vector / matrix of observed events whose elements should be named.
#' @param ec Vector / matrix of central exposure. The central exposure
#'   corresponds to the sum of the exposure period over the insured population.
#'   An individual experiencing an event of interest during the year will no
#'   longer be exposed afterwards and the exposure should be reduced
#'   accordingly.
#' @param lambda Smoothing parameter. If missing, an optimization procedure will
#'   be used to find the optimal smoothing parameter.
#' @param q Order of penalization. Polynoms of degrees `q - 1` are considered
#'   smooth and therefore unpenalized. The default of `2` should be suitable for
#'   most practical applications. Higher orders may cause numerical issues.
#' @param criterion Criterion to be used for the selection of the optimal
#'   smoothing parameter. Default is `"REML"` which stands for restricted
#'   maximum likelihood. Other options include `"AIC"`, `"BIC"` and `"GCV"`.
#' @param reg Should an approximate regression framework be used ?
#'   framework.
#' @param y Optional vector of observations whose elements should be named. Used
#'   only in the regression framework and even in this case will be
#'   automatically computed from the `d` and `ec` arguments if those are
#'   supplied. May be useful when using Whittaker-Henderson smoothing outside of
#'   the survival analysis framework.
#' @param wt Optional vector / matrix of weights. As for the observation vector
#'   / matrix `y`, used only in the regression framework and even in this case
#'   will be automatically computed if the `d` argument is supplied. May be
#'   useful when using Whittaker-Henderson smoothing outside of the survival
#'   analysis framework.
#' @param verbose Integer between 0 and 3. Control the level of informations
#'   that will be printed on screen during fitting.
#' @param ... Additional parameters passed to the smoothing function called.
#'
#' @returns An object of class `WH_1d` i.e. a list containing, among other
#'   things :
#' * `y` The observation vector/matrix, either supplied or computed as y = log(d) - log(ec)
#' * `y_hat` The vector/matrix of fitted value
#' * `std_y_hat` The vector/matrix of standard deviation associated with the fitted value
#' * `res` The vector/matrix of model deviance residuals
#' * `edf` The vector/matrix of effective degrees of freedom associated with each observation
#' * `diagnosis` A data.frame with one row containing the effective degrees of freedom of the model,
#'   the deviance of the fit as well as the AIC, BIC, GCV and REML criteria
#'
#' @examples
#'
#' d <- portfolio_mort$d
#' ec <- portfolio_mort$ec
#'
#' y <- log(d / ec)
#' y[d == 0 | ec == 0] <- NA
#' wt <- d
#'
#' # Maximum likelihood
#' WH(d, ec) # automatic smoothing parameter selection via REML
#' WH(d, ec, lambda = 1e2) # fixed smoothing parameter
#' WH(d, ec, criterion = "GCV") # alternative criterion for smoothing parameter selection
#'
#' # Regression
#' WH(y = y, wt = wt) # regression framework is default when y is supplied
#' WH(d, ec, reg = TRUE, lambda = 1e2) # forces computation of y from d and ec
#'
#' @export
WH <- function(d, ec, lambda = NULL, q = 2, criterion, reg, y, wt, verbose = 1, ...) {

  if (missing(d) && missing(y)) stop("An observation vector (either d or y) must be supplied")
  if (!missing(d) && !missing(y)) stop("Only one observation vector (either d or y) must be supplied")

  # Guess dimension of estimation problem
  target <- if (missing(d)) y else d
  dims <- if (is.null(dim(target))) 1 else 2
  if (dims > 2) stop("At the moment, WH only handles smoothing in 1 or 2 dimensions, sorry !")

  if (missing(reg)) reg <- !missing(y)

  # Check availability of observation and weight vectors
  if (reg) {
    if (missing(y)) {
      if (!missing(d) && !missing(ec)) {

        if (dims == 1) {

          if (!is.numeric(d)) stop("d should be a numeric vector")
          if (!is.numeric(ec)) stop("ec should be a numeric vector")
          if (length(d) != length(ec)) stop("Length of d and ec must match")

        } else {

          if (is.null(dim(d)) || length(dim(d)) != dims) stop("d should be a numeric matrix")
          if (is.null(dim(ec)) || length(dim(ec)) != dims) stop("ec should be a numeric matrix")
          if (any(dim(d) != dim(ec))) stop("Dimensions of d and ec must match")

        }
        if (verbose > 0) message("Computing y as log(d / ec)")
        y <- log(d / ec)
        y[d == 0 | ec == 0] <- NA
      } else {
        stop("Either y or both d and ec required in the regression framework")
      }
    } else {
      if (dims == 1) {
        if (!is.numeric(y)) stop("y should be a numeric vector")
      } else {
        if (is.null(dim(y)) || length(dim(y)) != 2) stop("y should be a numeric matrix")
      }
    }
    if (missing(wt)) {
      if (!missing(d)) {
        if (verbose > 0) message("Using d as weights")
        wt <- d
      } else {
        stop("Either wt or d needs to be supplied in the regression framework")
      }
    } else {
      if (dims == 1) {
        if (!is.numeric(wt)) stop("wt should be a numeric vector")
        if (length(y) != length(wt)) stop("Length of y and wt must match")
      } else {
        if (is.null(dim(wt)) || length(dim(wt)) != 2) stop("wt should be a numeric matrix")
        if (any(dim(y) != dim(wt))) stop("Dimensions of y and wt must match")
      }
    }
    if (dims == 1) {
      if (is.null(names(y))) {
        y_names <- seq_along(y)
        if (verbose > 0) warning("No names found for y, setting names to: ", paste(range(y_names), collapse = " - "))
      }
    } else {
      if (is.null(rownames(y))) {
        rownames(y) <- seq_along(dim(y)[[1]])
        if (verbose > 0) warning("No names found for y rows, setting names to: ", paste(range(rownames(y)), collapse = " - "))
      }
      if (is.null(colnames(y))) {
        colnames(y) <- seq_along(dim(y)[[2]])
        if (verbose > 0) warning("No names found for y columns, setting names to: ", paste(range(colnames(y)), collapse = " - "))
      }
    }
  } else {
    if (missing(d)) stop("d required in the maximum likelihood framework")
    if (missing(ec)) stop("ec required in the maximum likelihood framework")
    if (dims == 1) {
      if (!is.numeric(d)) stop("d should be a numeric vector")
      if (!is.numeric(ec)) stop("ec should be a numeric vector")
      if (length(d) != length(ec)) stop("Dimensions of d and ec must match")
      if (is.null(names(d))) {
        d_names <- seq_along(d)
        if (verbose > 0) warning("No names found for d, setting names to: ", paste(range(d_names), collapse = " - "))
      }
    } else {
      if (is.null(dim(d)) || length(dim(d)) != dims) stop("d should be a numeric matrix")
      if (is.null(dim(ec)) || length(dim(ec)) != dims) stop("ec should be a numeric matrix")
      if (any(dim(d) != dim(ec))) stop("Dimensions of d and ec must match")
      if (is.null(rownames(d))) {
        rownames(d) <- seq_along(dim(d)[[1]])
        if (verbose > 0) warning("No names found for d rows, setting names to: ", paste(range(rownames(d)), collapse = " - "))
      }
      if (is.null(colnames(d))) {
        colnames(d) <- seq_along(dim(d)[[2]])
        if (verbose > 0) warning("No names found for d columns, setting names to: ", paste(range(colnames(d)), collapse = " - "))
      }
    }
  }

  # Check method, criterion and smoothing parameters
  if (!is.null(lambda)) {
    if (length(lambda) > dims) stop(paste0("lambda should be of length ", dims))
    if (!is.numeric(lambda) || any(lambda < 0)) stop(
      "lambda should only have positive elements")
  }

  valid_criteria <- c("REML", "AIC", "BIC", "GCV")
  if (missing(criterion)) criterion <- "REML"
  if (length(criterion) != 1 || !(criterion %in% valid_criteria)) stop(
    "criterion should be one of: ", paste(valid_criteria, collapse = ", "))

  method <- if (is.null(lambda)) "outer" else "fixed_lambda"

  # Check penalization order
  if (!is.numeric(q) || length(q) > dims || length(q) == 0 || any(q <= 0) || max(abs(q - round(q))) > .Machine$double.eps ^ 0.5) stop(
    paste0("q should be an integer vector of length ", dims))
  if (any(q > 3) && verbose > 0) warning(
    "Differences of order q > 3 are not recommanded and may cause numerical instability")

  # Use the adequate function and arguments
  what <- paste("WH", method, sep = "_")

  args <- c(
    if (reg) list(y = y, wt = wt) else list(d = d, ec = ec),
    list(reg = reg, lambda = lambda, q = q, verbose = verbose),
    if (is.null(lambda)) list(criterion = criterion),
    list(...))

  do.call(what, args)
}

# Extrapolation----

#' Predict new values using a fitted 1D WH model
#'
#' Extrapolate the model for new observations.
#'
#' @param object An object of class `"WH_1d"` returned by the [WH()] function
#' @param newdata A vector containing the position of new observations.
#'   Observations from the fit will automatically be added to this, in the
#'   adequate order
#' @param ... Not used
#'
#' @returns An object of class `"WH_1d"` with additional components for model
#'   prediction.
#'
#' @examples
#'
#' object <- WH(portfolio_mort$d, portfolio_mort$ec)
#' object_extra <- predict(object, newdata = 40:99)
#' plot(object_extra)
#'
#' @export
predict.WH_1d <- function(object, newdata = NULL, ...) {

  if (!inherits(object, "WH_1d")) stop("object must be of class WH_1d")
  if (!is.numeric(newdata)) stop("newdata should be a vector containing the names of the predicted values")

  # Requires the following elements from object: lambda, q, wt, tXWz
  data <- as.numeric(names(object$wt))

  n <- length(data)
  n_tot <-length(newdata)
  ind_fit <- sum(newdata < min(data)) + seq_len(n)

  P <- create_P_compact_cpp(n_tot, object$q)
  lambda_P <- combine_lambda_P(object$lambda, P)

  wt <- tXWz <- rep(0, n_tot)
  wt[ind_fit] <- object$wt
  tXWz[ind_fit] <- object$tXWz

  fit <- fit_pen_reg(wt, tXWz, lambda_P)
  diag_V_pred <- diag_V_compact_cpp(fit$R)

  y_pred <- fit$y_hat
  std_y_pred <- sqrt(diag_V_pred)

  names(y_pred) <- names(std_y_pred) <- newdata

  object$y_pred <- y_pred
  object$std_y_pred <- std_y_pred
  object$R_pred <- fit$R

  return(object)
}

#' Predict new values using a fitted 2D WH model
#'
#' Extrapolate the model for new observations in a way that is consistent with
#' the fitted values
#'
#' @param object An object of class `"WH_2d"` returned by the [WH()] function
#' @param newdata A list containing two vectors indicating the new observation
#'   positions
#' @param ... Not used
#'
#' @returns An object of class `"WH_2d"` with additional components for model
#'   prediction.
#'
#' @examples
#' object <- WH(portfolio_LTC$d, portfolio_LTC$ec)
#' object_extra <- predict(object, newdata = list(age = 60:109, duration = 0:19))
#' plot(object_extra)
#'
#' @export
predict.WH_2d <- function(object, newdata = NULL, ...) {

  # Requires the following elements from object: lambda, q, wt, tXWz, R, y_hat, std_y_hat

  if (!inherits(object, "WH_2d")) stop("object must be of class WH_2d")
  if (length(newdata) != 2 || !is.numeric(newdata[[1]]) || !is.numeric(newdata[[2]])) stop(
    "newdata should be a list with two elements containing the row names and column names for predicted values")

  data <- dimnames(object$wt) |> lapply(as.numeric)

  n <- sapply(data, length)
  n_tot <- sapply(newdata, length)
  n_inf <- mapply(\(x,y) sum(y < min(x)), data, newdata)

  ind_mat <- matrix(seq_len(prod(n_tot)), n_tot[[1]])
  ind_rows <- list(n_inf[[1]] + seq_len(n[[1]]), n_inf[[2]] + seq_len(n[[2]]))
  ind_fit <- ind_mat[ind_rows[[1]], ind_rows[[2]]]; dim(ind_fit) <- NULL

  P <- mapply(create_P_compact_cpp, n_tot, object$q)
  lambda_P <- list(object$lambda[[1]] * P[[1]], object$lambda[[2]] * P[[2]])
  P_pred <- combine_P_compact(lambda_P)
  R_22 <- submatrix_compact_cpp(P_pred, ind_fit) |> cholesky_compact_lapack(inplace = TRUE)

  y_new <- get_y_new_compact(y = object$y_hat, lambda_P, R_22, n_tot, ind_fit)
  diag_V_pred <- get_diag_V_pred_compact(object$R, lambda_P, R_22, n, n_tot, ind_fit) + diag_V_compact_cpp(R_22)

  y_pred <- std_y_pred <- numeric(prod(n_tot))
  y_pred[ind_fit] <- object$y_hat
  y_pred[- ind_fit] <- y_new
  std_y_pred[ind_fit] <- object$std_y_hat
  std_y_pred[- ind_fit] <- sqrt(diag_V_pred)

  dim(y_pred) <- dim(std_y_pred) <- n_tot # set dimension for output matrices
  dimnames(y_pred) <- dimnames(std_y_pred) <- newdata # set names for output matrices

  object$y_pred <- y_pred
  object$std_y_pred <- std_y_pred
  object$lambda_P <- lambda_P
  object$R_22 <- R_22
  object$ind_fit <- ind_fit

  return(object)
}

# Variance-covariance matrix----

#' Compute variance-covariance matrix of fitted 1D WH model
#'
#' The variance-covariance matrix may be useful in case confidence intervals are
#' required for quantities derived from the fitted values.
#'
#' @param object An object of class `"WH_1d"` returned by the [WH()] function
#' @param pred Should the variance-covariance matrix include the extrapolated
#'   values as well (if any) ?
#' @param ... Not used
#'
#' @returns The variance-covariance matrix for the fitted values
#'
#' @examples
#'
#' object <- WH(portfolio_mort$d, portfolio_mort$ec)
#' vcov(object)
#'
#' object_extra <- predict(object, newdata = 40:99)
#' V <- vcov(object_extra)
#'
#' @export
vcov.WH_1d <- function(object, pred = TRUE, ...) {

  pred <- pred && !is.null(object$R_pred)
  R <- if (pred) object$R_pred else object$R
  invert_cholesky_compact_lapack(R) |> tcrossprod()
}

#' Compute variance-covariance matrix of fitted 1D WH model
#'
#' The variance-covariance matrix may be useful in case confidence intervals are
#' required for quantities derived from the fitted values.
#'
#' @param object An object of class `"WH_2d"` returned by the [WH()] function
#' @param pred Should the variance-covariance matrix include the extrapolated
#'   values as well (if any) ?
#' @param ... Not used
#'
#' @returns The variance-covariance matrix for the fitted values
#'
#' @examples
#'
#' object <- WH(portfolio_LTC$d, portfolio_LTC$ec)
#' V <- vcov(object)
#'
#' object_extra <- predict(object, newdata = list(age = 60:109, duration = 0:19))
#' V <- vcov(object_extra)
#'
#' @export
vcov.WH_2d <- function(object, pred = TRUE, ...) {

  pred <- pred && !is.null(object$R_22)
  if (pred) {

    prod_n_tot <- prod(dim(object$y_pred))

    K <- invert_cholesky_compact_lapack(object$R)
    V_11 <- tcrossprod(K)

    K_aug <- numeric(prod_n_tot * prod_n_tot); dim(K_aug) <- rep(prod_n_tot, 2)
    K_aug[object$ind_fit, object$ind_fit] <- K

    P_aux <- get_prod_P_K_compact_cpp(K_aug, object$lambda_P)[- object$ind_fit, object$ind_fit]
    A <- backsolve_mat_compact_lapack(object$R_22, backsolve_mat_compact_lapack(object$R_22, P_aux, transpose = TRUE))
    K_22 <- invert_cholesky_compact_lapack(object$R_22)

    V_22 <- tcrossprod(A) + tcrossprod(K_22)
    V_12 <- - backsolve_mat_compact_lapack(object$R, t(A))

    V_pred <- K_aug
    V_pred[object$ind_fit, object$ind_fit] <- V_11
    V_pred[object$ind_fit, - object$ind_fit] <- V_12
    V_pred[- object$ind_fit, object$ind_fit] <- t(V_12)
    V_pred[- object$ind_fit, - object$ind_fit] <- V_22
    V_pred

  } else {
    invert_cholesky_compact_lapack(object$R) |> tcrossprod()
  }
}

# Formatting----

#' Store WH model fit results in a data.frame
#'
#' @param object An object of class  `"WH_1d"` or `"WH_2d"` returned
#'   by the [WH()] function
#' @param dim1 The (optional) name to be given to the first dimension
#' @param dim2 The (optional) name to be given to the second dimension
#'
#' @returns A data.frame gathering information about the fitted and predicted
#'   values, the model variance, residuals and effective degrees of freedom...
#'
#' @export
output_to_df <- function(object, dim1 = "x", dim2 = "z") {

  if (!inherits(object, c("WH_1d", "WH_2d"))) stop("object must be of class WH_1d or WH_2d")
  if (length(dim1) != 1 || length(dim2) != 1) stop("The dim1 and dim2 optional arguments should be of length 1")

  if ("y_pred" %in% names(object)) {
    object$y_hat <- object$y_pred
    object$std_y_hat = object$std_y_pred
  }

  df <- data.frame(y_hat = c(object$y_hat), std_y_hat = c(object$std_y_hat))

  if (inherits(object, "WH_1d")) {

    x <- as.numeric(names(object$y))
    x_pred <- as.numeric(names(object$y_hat))

    df[[dim1]] <- x_pred

    df$y <- df$wt <- df$res <- df$edf <- NA_real_
    df$y[x_pred %in% x] <- c(object$y)
    df$wt[x_pred %in% x] <- c(object$wt)
    df$res[x_pred %in% x] <- c(object$res)
    df$edf[x_pred %in% x] <- c(object$edf)

    if("ec" %in% names(object) && "d" %in% names(object)) {

      df$ec <- df$d <- NA_real_
      df$ec[x_pred %in% x] <- c(object$ec)
      df$d[x_pred %in% x] <- c(object$d)
    }

    df[, c(dim1, intersect(names(object), c("y", "y_hat", "std_y_hat", "res", "edf")))]

  } else {

    x <- as.numeric(rownames(object$y))
    x_pred <- as.numeric(rownames(object$y_hat))

    z <- as.numeric(colnames(object$y))
    z_pred <- as.numeric(colnames(object$y_hat))

    df[[dim1]] <- rep(x_pred, times = length(z_pred))
    df[[dim2]] <- rep(z_pred, each = length(x_pred))

    df$y <- df$wt <- df$res <- df$edf <- NA_real_
    df$y[df[[dim1]] %in% x & df[[dim2]] %in% z] <- c(object$y)
    df$wt[df[[dim1]] %in% x & df[[dim2]] %in% z] <- c(object$wt)
    df$res[df[[dim1]] %in% x & df[[dim2]] %in% z] <- c(object$res)
    df$edf[df[[dim1]] %in% x & df[[dim2]] %in% z] <- c(object$edf)

    if("ec" %in% names(object) && "d" %in% names(object)) {

      df$ec <- df$d <- NA_real_
      df$ec[df[[dim1]] %in% x & df[[dim2]] %in% z] <- c(object$ec)
      df$d[df[[dim1]] %in% x & df[[dim2]] %in% z] <- c(object$d)
    }

    df[, c(dim1, dim2, intersect(names(object), c("y", "y_hat", "std_y_hat", "res", "edf")))]
  }
}

# Display----

#' Display of 1D WH object
#'
#' @param x An object of class `"WH_1d"` returned by the [WH()] function
#' @param ... Not used
#'
#' @returns Invisibly returns `x`.
#'
#' @examples
#' WH(portfolio_mort$d, portfolio_mort$ec)
#'
#' @export
print.WH_1d <- function(x, ...) {

  if (!inherits(x, "WH_1d")) stop("x must be of class WH_1d")

  cat("An object fitted using the WH function\n")
  cat("Initial data contains", length(x$y), "data points:\n")
  cat("  Observation positions: ", as.numeric(names(x$y)[[1]]), " to ",as.numeric(names(x$y)[[length(x$y)]]), "\n")
  cat("Smoothing parameter selected:", format(x$lambda, digits = 2), "\n")
  cat("Associated degrees of freedom:", format(sum(x$edf), digits = 2), "\n\n")
  invisible(x)
}

#' Display of 2D WH object
#'
#' @param x An object of class `"WH_2d"` returned by the [WH()] function
#' @param ... Not used
#'
#' @returns Invisibly returns `x`.
#'
#' @examples
#' WH(portfolio_LTC$d, portfolio_LTC$ec)
#'
#' @export
print.WH_2d <- function(x, ...) {

  if (!inherits(x, "WH_2d")) stop("x must be of class WH_2d")

  cat("An object fitted using the WH function\n")
  cat("Initial data contains", prod(dim(x$y)), "data points:\n")
  cat("  First dimension: ", as.numeric(rownames(x$y)[[1]]), " to ",as.numeric(rownames(x$y)[[nrow(x$y)]]), "\n")
  cat("  Second dimension: ", as.numeric(colnames(x$y)[[1]]), " to ",as.numeric(colnames(x$y)[[ncol(x$y)]]), "\n")
  cat("Smoothing parameters selected:", format(x$lambda, digits = 2), "\n")
  cat("Associated degrees of freedom:", format(sum(x$edf), digits = 2), "\n\n")
  invisible(x)
}

# Plots----

#' Plot 1D WH fit
#'
#' @param x An object of class `"WH_1d"` returned by the [WH()] function
#' @param what What should be plotted. Should be one of `fit` (the
#'   default), `res` for residuals and `edf` for the effective degrees
#'   of freedom.
#' @param trans An (optional) transformation to be applied to the data. By
#'   default the identity function
#' @param ... Not used
#'
#' @returns A plot representing the desired element of the fit
#'
#' @examples
#' d <- portfolio_mort$d
#' ec <- portfolio_mort$ec
#'
#' WH(d, ec) |> plot()
#' WH(d, ec) |> plot("res")
#' WH(d, ec) |> plot("edf")
#'
#' @export
plot.WH_1d <- function(x, what = "fit", trans, ...) {

  if (!inherits(x, "WH_1d")) stop("x must be of class WH_1d")
  whats <- c("fit", "res", "edf")
  if (length(what) != 1 || !(what %in% whats)) stop(
    "what should be one of ", paste(whats, collapse = ", "))

  df <- output_to_df(x)
  if (missing(trans)) trans <- \(x) x

  switch(
    what,
    fit = {
      plot(df$x, trans(df$y),
           xlab = "x",
           ylab = "log - hazard rate",
           xlim = range(df$x),
           ylim = trans(range(c(df$y_hat - 2 * df$std_y_hat),
                              c(df$y_hat + 2 * df$std_y_hat))))
      graphics::lines(df$x, trans(df$y_hat), col = "blue")
      graphics::lines(df$x, trans(df$y_hat - 2 * df$std_y_hat), col = "red", lty = 3)
      graphics::lines(df$x, trans(df$y_hat + 2 * df$std_y_hat), col = "red", lty = 3)
    },
    res = {
      plot(df$x, df$res, xlab = "x", ylab = "deviance residuals", type = "b")
      graphics::abline(a = 0, b = 0, lty = 2, col = "blue")
    },
    edf = {
      plot(df$x, df$edf, xlab = "x", ylab = "degrees of freedom", type = "b")
      graphics::abline(a = 0, b = 0, lty = 2, col = "blue")
      graphics::abline(a = 1, b = 0, lty = 2, col = "blue")
    })
}

#' Plot 2D WH fit
#'
#' @param x An object of class `"WH_2d"` returned by the [WH()] function
#' @param what What should be plotted (y_hat, std_y_hat, res, edf)
#' @param trans An (optional) transformation to be applied to the data
#' @param ... Not used
#'
#' @returns A plot representing the given element of the fit...
#'
#' @examples
#' d  <- portfolio_LTC$d
#' ec <- portfolio_LTC$ec
#'
#' WH(d, ec) |> plot()
#' WH(d, ec) |> plot("std_y_hat")
#'
#' @export
plot.WH_2d <- function(x, what = "y_hat", trans, ...) {

  if (!inherits(x, "WH_2d")) stop("x must be of class WH_2d")
  whats <- c("y_hat", "std_y_hat", "res", "edf")
  if (length(what) != 1 || !(what %in% whats)) stop(
    "what should be one of ", paste(whats, collapse = ", "))

  df <- output_to_df(x)
  if (missing(trans)) trans <- \(x) x

  x <- unique(df$x)
  z <- unique(df$z)

  data <- matrix(c(df[[what]]), length(z), length(x), byrow = TRUE)

  graphics::contour(
    z, x, data,
    nlevels = 20,
    xlab = "x",
    ylab = "z",
    main = switch(
      what,
      y_hat = "log - hazard rate",
      std_y_hat = " standard deviation of log - hazard rate",
      edf = "degrees of freedom"))
}

# Fit functions----

fit_pen_reg <- function(wt, tXWz, lambda_P) {

  lambda_P[nrow(lambda_P), ] <- lambda_P[nrow(lambda_P), ] + wt
  R <- cholesky_compact_lapack(lambda_P, inplace = TRUE)
  y_hat <- backsolve2(R, tXWz) # O(p ^ 2) cost, much faster than inverting Psi
  dim(y_hat) <- dim(wt)

  list(y_hat = y_hat, R = R)
}

fit_WH_fixed_lambda_reg <- function(y, wt, tXWz, lambda_P, q, lambda) {

  fit <- fit_pen_reg(wt, tXWz, lambda_P)

  res <- compute_res_reg(y, fit$y_hat, wt)
  dev <- sum(res * res)

  RESS <- get_RESS(fit$y_hat, q)
  pen <- sum(lambda * RESS)

  list(wt = wt, tXWz = tXWz, y_hat = fit$y_hat, R = fit$R, res = res, dev = dev, RESS = RESS, pen = pen)
}

get_valid_beta <- function(y_hat, new_y_hat, d, off, q, lambda, counter, verbose, tol_dev, old_dev_pen, step_halving, halving_tr = 1024) {

  if (counter > 1) step_y <- new_y_hat - y_hat

  step_factor <- 1
  valid_candidate <- FALSE

  while (!valid_candidate) {

    y_hat_candidate <- if (counter > 1 && step_factor < halving_tr) (y_hat + step_y / step_factor) else y_hat

    new_wt <- exp(y_hat_candidate + off)
    res <- compute_res_deviance(d, new_wt)
    dev <- sum(res * res)

    RESS <- get_RESS(y_hat_candidate, q)
    pen <- sum(lambda * RESS)

    dev_pen <- dev + pen

    valid_candidate <- counter == 1 || !step_halving || dev_pen <= old_dev_pen

    if (verbose > 1) print_fit_infos(
      "PIRLS", lambda, counter, "penalized deviance", dev_pen, old_dev_pen,
      if (counter != 1 && step_halving) valid_candidate, step_factor)

    if (!valid_candidate) step_factor <- 2 * step_factor
  }
  if (verbose > 0) {
    if (step_factor >= halving_tr) message(paste0(
      "PIRLS step divider reached ", halving_tr, ". Reverting to initial parameter estimate"))
    if (!valid_candidate) warning("PIRLS failed to improve penalized deviance")
  }
  list(y_hat = y_hat_candidate, new_wt = new_wt, res = res, dev = dev, RESS = RESS, pen = pen, dev_pen = dev_pen)
}

fit_WH_fixed_lambda_ml <- function(d, off, lambda_P, q, lambda, verbose, tol_dev, eta_start, step_halving) {

  # Use the saturated model rates as starting estimates and compute associated weights
  y_hat <- if (is.null(eta_start)) init_y_hat(d, off) else eta_start  # the small constant ensures this is always defined
  new_wt <- exp(y_hat + off)

  counter <- 0
  dev_pen <- Inf
  cond <- TRUE

  while (cond) {

    counter <- counter + 1
    old_dev_pen <- dev_pen

    wt <- new_wt
    z <- compute_z(y_hat, new_wt, d)
    tXWz <- wt * z; dim(tXWz) <- NULL

    fit <- fit_pen_reg(wt = wt, tXWz = tXWz, lambda_P)

    if (counter == 1) y_hat <- fit$y_hat

    valid_beta <- get_valid_beta(
      y_hat, new_y_hat = fit$y_hat, d, off, q, lambda, counter, verbose, tol_dev, old_dev_pen, step_halving)

    y_hat <- valid_beta$y_hat
    new_wt <- valid_beta$new_wt
    dev_pen <- valid_beta$dev_pen

    cond <- old_dev_pen - dev_pen >= tol_dev * old_dev_pen
  }
  if (verbose > 0) print_fit_infos("PIRLS", lambda, counter, "penalized deviance", dev_pen, old_dev_pen, final = TRUE)

  list(wt = wt, tXWz = tXWz, z = z, y_hat = y_hat, R = fit$R,
       res = valid_beta$res, dev = valid_beta$dev, RESS = valid_beta$RESS, pen = valid_beta$pen)
}

# Main functions----

WH_fixed_lambda <- function(
    d, ec, y, wt, reg = missing(d), lambda = 1e3, q = 2,
    verbose = 1, tol_dev = .Machine$double.eps^.5, eta_start = NULL, step_halving = TRUE) {

  # Guess dimension of estimation problem
  target <- if (reg) y else d
  dims <- if (is.null(dim(target))) 1 else length(dim(target))
  n <- if (dims == 1) length(target) else dim(target)

  # Extend arguments to multidimensional
  lambda <- extend_vector(lambda, dims)
  q <- extend_vector(q, dims)

  # Retrieve model matrix and diagonal of penalization (diagonal) matrix
  model_matrices <- get_model_matrices(n, q)
  P <- model_matrices$P
  s <- model_matrices$s

  lambda_P <- combine_lambda_P(lambda, P)

  which_pos <- if (reg) which(wt != 0) else which(ec != 0) # indices of observation with non-zero weight
  n_pos <- length(which_pos) # number of observations with non-zero weight

  if (reg) {

    fit <- fit_WH_fixed_lambda_reg(
      y = y, wt = wt, tXWz = wt * y, lambda_P = lambda_P, q = q, lambda = lambda)
    z <- y

  } else {

    off <- ifelse(ec == 0, - Inf, log(ec))
    y <- ifelse(d == 0 | ec == 0, NA, log(d) - off)

    fit <- fit_WH_fixed_lambda_ml(
      d = d, off = off, lambda_P = lambda_P, q = q, lambda = lambda,
      verbose = verbose, tol_dev = tol_dev, step_halving = step_halving, eta_start = eta_start)
    wt <- fit$wt
    z <- fit$z
  }

  res <- fit$res

  diag_V <- diag_V_compact_cpp(fit$R)
  l_sat <- if (reg) get_l_sat(n = n_pos) else 0

  y_hat <- fit$y_hat # fit results (including predicted values for zero-weight observations)
  std_y_hat <- sqrt(diag_V) # standard deviation of the fit

  edf <- wt * diag_V # effective degrees of freedom by observation
  sum_edf <- sum(edf) # effective degrees of freedom

  diagnosis <- get_diagnosis(
    dev = fit$dev, pen = fit$pen, sum_edf = sum_edf, n_pos = n_pos,
    p = n, q = q, s = s, lambda = lambda, R = fit$R, l_sat = l_sat)

  # Output vector / matrix names
  if (dims == 1) {

    names(y_hat) <- names(std_y_hat) <- names(res) <- names(edf) <-
      names(wt) <- names(z) <- names(y) # set names for output vectors

  } else {

    dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(edf) <-
      dim(wt) <- dim(z) <- dim(y) # set dimensions for output matrices
    dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(edf) <-
      dimnames(wt) <- dimnames(z) <- dimnames(y) # set names for output matrices
  }

  out <- list(
    y = y, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf = edf, diagnosis = diagnosis,
    wt = wt, tXWz = fit$tXWz, R = fit$R, lambda = lambda, q = q)
  class(out) <- if (dims == 1) "WH_1d" else "WH_2d"

  return(out)
}

WH_outer <- function(
    d, ec, y, wt, reg = missing(d), q = 2, criterion = "REML",
    lambda = NULL, verbose = 1, tol_crit = .Machine$double.eps^.5, tol_dev = .Machine$double.eps^.5,
    step_halving = c(TRUE, TRUE), lambda_start_ratio = 0.2, auto_transpose = TRUE) {

  # Guess dimension of estimation problem
  target <- if (reg) y else d
  dims <- if (is.null(dim(target))) 1 else length(dim(target))
  n <- if (dims == 1) length(target) else dim(target)

  # Extend arguments to multidimensional
  lambda <- extend_vector(lambda, dims)
  q <- extend_vector(q, dims)

  auto_transpose <- auto_transpose && (dims > 1) && ((q[[2]] + 1) * n[[1]] > (q[[1]] + 1) * n[[2]])
  if (auto_transpose) {
    if (reg) {
      y <- t(y)
      wt <- t(wt)
    } else {
      d <- t(d)
      ec <- t(ec)
    }
    n <- rev(n)
    lambda <- rev(lambda)
    q <- rev(q)
  }

  criterion_display <- if (criterion == "REML" & !reg) "LAML" else criterion
  eta_start <- NULL

  # Retrieve model matrix and diagonal of penalization (diagonal) matrix
  model_matrices <- get_model_matrices(n, q)
  P <- model_matrices$P
  s <- model_matrices$s

  which_pos <- if (reg) which(wt != 0) else which(ec != 0) # indices of observation with non-zero weights
  n_pos <- length(which_pos) # number of observations with non-zero weight

  if (reg) {
    tXWz <- wt * y
  } else {
    off <- ifelse(ec == 0, - Inf, log(ec))
  }

  if (is.null(lambda) && dims == 2) {

    diag_P <- if (dims == 1) {
      P[nrow(P), ]
    } else {
      diag_P1 <- P[[1]][nrow(P[[1]]), ]
      diag_P2 <- P[[2]][nrow(P[[2]]), ]
      list(rep(diag_P1, n[[2]]), rep(diag_P2, each = n[[1]]))
    }
    if (!reg) wt <- exp(init_y_hat(d, off) + off)
    lambda <- find_starting_lambda(wt = wt, diag_P = diag_P, lambda_start_ratio)
  }

  counter_outer <- 0
  old_score <- Inf
  cond <- TRUE
  lambda_track <- list()

  # Auxiliary function for smoothing parameter optimization
  WH_outer_aux <- function(log_lambda) {

    counter_outer <<- counter_outer + 1 # update iteration counter (globally)
    lambda <- exp(log_lambda)
    lambda_track[[counter_outer]] <<- lambda

    lambda_P <- combine_lambda_P(lambda, P)

    if (reg) {

      fit <- fit_WH_fixed_lambda_reg(
        y = y, wt = wt, tXWz = tXWz, lambda_P = lambda_P, q = q, lambda = lambda)

    } else {

      fit <- fit_WH_fixed_lambda_ml(
        d = d, off = off, lambda_P = lambda_P, q = q, lambda = lambda,
        verbose = verbose - 1, tol_dev = tol_dev, eta_start = eta_start, step_halving = step_halving[[2]])
      eta_start <<- fit$y_hat  # update starting parameters (globally)
    }
    l_sat <- if (reg) get_l_sat(n = n_pos) else 0

    if (criterion != "REML") sum_edf <- get_edf_tot(R = fit$R, wt = fit$wt) # effective degrees of freedom
    score <- switch(
      criterion,
      AIC = fit$dev + 2 * sum_edf,
      BIC = fit$dev + log(n_pos) * sum_edf,
      GCV = n_pos * fit$dev / (n_pos - sum_edf) ^ 2,
      REML = compute_REML(dev = fit$dev, pen = fit$pen, p = n, q = q, s = s, lambda = lambda, R = fit$R, l_sat = l_sat)
    )
    if (verbose > 1) print_fit_infos("optim", lambda, counter_outer, criterion_display, score, old_score)
    old_score <<- score # update score (globally)
    return(score)
  }

  # Call the optimization process to find the optimal smoothing parameter
  lambda <- if (dims == 1) {
    exp(stats::optimize(f = WH_outer_aux, interval = c(- 6, 12) * log(10))$minimum) # Brent (1973) algorithm
  } else {
    exp(stats::optim(par = log(lambda), fn = WH_outer_aux, control = list(reltol = tol_crit))$par) # Nelder and Mead (1965) algorithm
  }

  if (auto_transpose) {
    if (reg) {
      y <- t(y)
      wt <- t(wt)
    } else {
      d <- t(d)
      ec <- t(ec)
      eta_start <- t(eta_start)
    }
    lambda <- rev(lambda)
    q <- rev(q)
  }

  # Call to a fixed smoothing parameter to get all model statistics
  out <- WH_fixed_lambda(
    d = d, ec = ec, y = y, wt = wt, reg = reg, lambda = lambda, q = q,
    verbose = verbose - 1, tol_dev = tol_dev, eta_start = eta_start, step_halving = step_halving[[1]])

  # Print summary of the optimization process
  if (verbose > 0) print_fit_infos("Outer", lambda, counter_outer, criterion_display, out$diagnosis[[criterion]], final = TRUE)

  out$lambda_track <- lambda_track
  return(out)
}


