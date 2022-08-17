# Wrappers functions----

#' 1D Whittaker-Henderson Smoothing
#'
#' Main package function to apply Whittaker-Henderson smoothing in a
#' unidimensional survival analysis framework. It takes as input a vector of
#' observed events and a vector of associated central exposure, both depending
#' on a single covariate, and build a smooth version of the log-hazard rate.
#' Smoothing parameters may be supplied or automatically chosen according to a
#' specific criterion such as \code{"REML"} (the default), \code{"AIC"}, \code{"BIC"} or
#' \code{"GCV"}. Whittaker-Henderson may be applied in a full maximum likelihood
#' framework or an asymptotic (approximate) gaussian framework.
#'
#' @param d Vector of observed events whose elements should be named.
#' @param ec Vector of central exposure. The central exposure corresponds to the
#'   sum of the exposure period over the insured population. An individual
#'   experiencing an event of interest during the year will no longer be exposed
#'   afterward and the exposure should be reduced accordingly.
#' @param lambda Smoothing parameter. If missing, an optimization procedure will
#'   be used to find the optimal smoothing parameter. If supplied, no optimal
#'   smoothing parameter search will take place unless the \code{criterion}
#'   argument is also supplied, in which case \code{lambda} will be used as the
#'   starting parameter for the optimization procedure.
#' @param criterion Criterion to be used for the selection of the optimal
#'   smoothing parameter. Default is \code{"REML"} which stands for restricted
#'   maximum likelihood. Other options include \code{"AIC"}, \code{"BIC"} and \code{"GCV"}.
#' @param method Method to be used to find the optimal smoothing parameter.
#'   Default to \code{"fixed_lambda"} if \code{lambda} is supplied, meaning no
#'   optimization is performed. Otherwise, if \code{criterion = "REML"} or the
#'   \code{criterion} argument is missing, default to \code{"fs"} which means the
#'   generalized Fellner-Schall method is used. For other criteria default to
#'   \code{"optim"} meaning the \code{optimize} function from package \code{stats}
#'   will be used.
#' @param q Order of penalization. Polynoms of degrees \code{q - 1} are
#'   considered smooth and are therefore unpenalized. Should be left to the
#'   default of \code{2} for most practical applications.
#' @param framework Default framework is \code{"ml"} which stands for maximum
#'   likelihood unless the \code{y} argument is also provided, in which case the
#'   \code{"reg"} or regression framework is used. The regression framework is an
#'   asymptotical approximation of the maximum likelihood likelihood framework.
#' @param y Optional vector of observations whose elements should be named. Used
#'   only in the regression framework and even in this case will be
#'   automatically computed from the \code{d} and \code{ec} arguments if those
#'   are supplied. May be useful when using Whittaker-Henderson smoothing
#'   outside of the survival analysis framework.
#' @param wt Optional vector of weights. As for the observation vector \code{y},
#'   used only in the regression framework and even in this case will be
#'   automatically computed if the \code{d} argument is supplied. May be useful
#'   when using Whittaker-Henderson smoothing outside of the survival analysis
#'   framework.
#' @param ... Additional parameters passed to the smoothing function called.
#'
#' @returns An object of class \code{WH_1d} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#'
#' @examples
#' keep <- which(portfolios_mort[[1]]$ec > 0)
#' d <- portfolios_mort[[1]]$d[keep]
#' ec <- portfolios_mort[[1]]$ec[keep]
#'
#' y <- log(d / ec)
#' y[d == 0] <- - 20
#' wt <- d
#'
#' # Maximum likelihood
#' WH_1d(d, ec, lambda = 1e2)
#' WH_1d(d, ec) # default generalized Fellner-Schall method
#' WH_1d(d, ec, method = "optim") # alternative method base on optimize function
#'
#' testthat::expect_equal(WH_1d(d, ec, method = "fs"),
#'                        WH_1d(d, ec, method = "optim"),
#'                        tolerance = 1e-1)
#' # generalized Fellner-Schall method is approximate in maximum likelihood framework
#'
#' WH_1d(d, ec, criterion = "GCV")
#' # alternative optimization criteria for smoothing parameter selection
#'
#' # Regression
#' WH_1d(y = y, wt = wt, lambda = 1e2)
#' # regression framework is default when y is supplied
#' WH_1d(d, y = y, lambda = 1e2) # d is used as wt
#' WH_1d(d, ec, framework = "reg", lambda = 1e2)
#' # Setting framework = "reg" forces computation of y from d and ec
#'
#' WH_1d(y = y, wt = wt)
#' WH_1d(y = y, wt = wt, method = "optim")
#' testthat::expect_equal(WH_1d(y = y, wt = wt, method = "fs"),
#'                        WH_1d(y = y, wt = wt, method = "optim"),
#'                        tolerance = 1e-6)
#' # generalized Fellner-Schall method is exact in regression framework
#'
#' WH_1d(y = y, wt = wt, criterion = "GCV")
#' # alternative optimization criteria for smoothing parameter selection
#'
#' @export
WH_1d <- function(d, ec, lambda, criterion, method, q = 2, framework, y, wt, ...) {

  if (missing(framework)) framework <- if (!missing(y)) "reg" else "ml"
  if (framework == "reg") {
    if (missing(y) && !missing(d) && !missing(ec)) {
      y <- log(d / ec)
      y[d == 0] <- - 20
    }
    if (missing(wt) && !missing(d)) wt <- d
  }

  if (missing(criterion)) criterion <- "REML"
  if (missing(method)) {
    method <- if (!missing(lambda)) "fixed_lambda" else {
      if (criterion != "REML") "optim" else "fs"
    }
  }
  if (missing(lambda)) lambda <- 1e3

  what <- paste("WH_1d", framework, method, sep = "_")
  args <- (if (framework == "reg") list(y = y, wt = wt) else list(d = d, ec = ec)) |>
    c(list(q = q, lambda = lambda),
      if (method == "optim") list(criterion = criterion),
      list(...))
  do.call(what, args)
}

#' 2D Whittaker-Henderson Smoothing
#'
#' #' Main package function to apply Whittaker-Henderson smoothing in a
#' bidimensional survival analysis framework. It takes as input a matrix of
#' observed events and a matrix of associated central exposure, both depending
#' on two covariates, and build a smooth version of the log-hazard rate.
#' Smoothing parameters may be supplied or automatically chosen according to a
#' specific criterion such as \code{"REML"} (the default), \code{"AIC"},
#' \code{"BIC"} or \code{"GCV"}. Whittaker-Henderson may be applied in a full
#' maximum likelihood framework or an asymptotic (approximate) gaussian
#' framework. As Whittaker-Henderson smoothing relies on full-rank smoothers,
#' computation time and memory usage in the bidimensional case may be
#' overwhelming and the function integrates an ad hoc rank-reduction procedure
#' to avoid such issues.
#'
#' @inheritParams WH_1d
#' @param d Matrix of observed events whose rows and columns should be named.
#'   Required in case of maximum likelihood estimation
#' @param ec Matrix of central exposure. Required in case of maximum likelihood
#'   estimation. The central exposure corresponds to the sum of the exposure
#'   period over the insured population. An individual experiencing an event of
#'   interest during the year will no longer be exposed afterward and the
#'   exposure should be reduced accordingly.
#' @param lambda Smoothing parameter vector of size \code{2}. If missing, an
#'   optimization procedure will be used to find the optimal smoothing
#'   parameter. If provided, no optimal smoothing parameter search will take
#'   place unless the \code{criterion} argument is also provided, in which case
#'   \code{lambda} will be used as the starting parameter for the optimization
#'   procedure.
#' @param p Optional vector of size \code{2}. Maximum number of eigenvectors to
#'   keep on each dimension after performing the eigen decomposition of the
#'   penalization matrix. If missing, will be automatically computed so that the
#'   dimensions of (square) matrices involved in the optimization problem
#'   remains lower that the \code{max_dim} argument
#' @param max_dim Maximal number of rows (or equivalently columns) of the square
#'   matrices that should be involved in the optimization problem. Default is
#'   \code{500}. Values higher than \code{2000} may cause issues with both
#'   memory usage and computation times.
#' @param q Order of penalization vector of size \code{2}. Polynoms of degrees
#'   \code{(q[[1]] - 1,q[[2]] - 1)} are considered smooth and are therefore
#'   unpenalized. Should be left to the default of \code{c(2,2)} for most
#'   practical applications.
#' @param y Optional matrix of observations whose rows and columns should be
#'   named. Used only in the regression framework and even in this case will be
#'   automatically computed if the \code{d} and \code{ec} arguments are
#'   supplied. May be useful when using Whittaker-Henderson smoothing outside of
#'   the survival analysis framework.
#' @param wt Optional matrix of weights. As for the observation vector \code{y},
#'   used only in the regression framework and even in this case will be
#'   automatically computed if the \code{d} argument is supplied. May be useful
#'   when using Whittaker-Henderson smoothing outside of the survival analysis
#'   framework.
#'
#' @returns An object of class \code{WH_2d} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#'
#' @examples
#' keep_age <- which(rowSums(portfolios_LTC[[1]]$ec) > 1e2)
#' keep_duration <- which(colSums(portfolios_LTC[[1]]$ec) > 1e2)
#'
#' d  <- portfolios_LTC[[1]]$d[keep_age, keep_duration]
#' ec <- portfolios_LTC[[1]]$ec[keep_age, keep_duration]
#'
#' y <- log(d / ec) # observation vector
#' y[d == 0] <- - 20
#' wt <- d
#'
#' # Maximum likelihood
#' WH_2d(d, ec, lambda = 1e2)
#' fit_fs <- WH_2d(d, ec) # default generalized Fellner-Schall method
#' fit_optim <- WH_2d(d, ec, method = "optim") # alternative method base on optim function
#' testthat::expect_equal(fit_fs, fit_optim, tolerance = 1e-1)
#' # generalized Fellner-Schall method is approximate in maximum likelihood framework
#'
#' WH_2d(d, ec, criterion = "GCV")
#' # alternative optimization criteria for smoothing parameter selectio
#'
#' # Regression
#' WH_2d(y = y, wt = wt, lambda = 1e2)
#' # regression framework is triggered when y is supplied
#' WH_2d(d, y = y, lambda = 1e2) # d is used as wt
#' WH_2d(d, ec, framework = "reg", lambda = 1e2)
#' # setting framework = "reg" forces computation of y from d and ec
#'
#' fit_fs <- WH_2d(y = y, wt = wt, method = "fs") # default generalized Fellner-Schall method
#' fit_optim <- WH_2d(y = y, wt = wt, method = "optim") # alternative method base on optim function
#' testthat::expect_equal(fit_fs, fit_optim, tolerance = 1e-4)
#' # generalized Fellner-Schall method is exact in regression framework
#'
#' WH_2d(y = y, wt = wt, criterion = "GCV")
#' # alternative optimization criteria for smoothing parameter selection
#'
#' keep_age2 <- which(rowSums(portfolios_LTC[[1]]$ec) > 0)
#' keep_duration2 <- which(colSums(portfolios_LTC[[1]]$ec) > 0)
#'
#' # Rank reduction
#' d  <- portfolios_LTC[[1]]$d[keep_age2, keep_duration2]
#' ec <- portfolios_LTC[[1]]$ec[keep_age2, keep_duration2]
#'
#' prod(dim(d)) # problem dimension is 1,262 !
#' WH_2d(d, ec)
#' # rank-reduction may be used to quickly find an approximate solution
#'
#' @export
WH_2d <- function(d, ec, lambda, criterion, method, p, max_dim = 250,
                  q = c(2, 2), framework, y, wt, ...) {

  if (missing(framework)) framework <- if (!missing(y)) "reg" else "ml"
  if (framework == "reg") {
    if (missing(y) && !missing(d) && !missing(ec)) {
      y <- log(d / ec)
      y[d == 0] <- - 20
    }
    if (missing(wt) && !missing(d)) wt <- d
  }

  if (missing(criterion)) criterion <- "REML"
  if (missing(method)) {
    method <- if (!missing(lambda)) "fixed_lambda" else {
      if (criterion != "REML") "optim" else "fs"
    }
  }
  if (missing(lambda)) lambda <- c(1e3, 1e3)

  dims <- if (!missing(y)) dim(y) else dim(d)
  if (missing(p)) {

    asp <- dims[[2]] / dims[[1]]
    max_ratio <- sqrt(max_dim / asp) / dims[[1]]
    p <- floor(pmin(max_ratio, 1) * dims)
  }

  what <- paste("WH_2d", framework, method, sep = "_")
  args <- (if (framework == "reg") list(y = y, wt = wt) else list(d = d, ec = ec)) |>
    c(list(p = p, q = q, lambda = lambda),
      if (method == "optim") list(criterion = criterion),
      list(...))
  do.call(what, args)
}

# Extrapolation----

#' Predict new values using a fitted 1D WH model
#'
#' Extrapolate the model for new values of the covariates
#'
#' @param object An object of class \code{"WH_1d"} returned by the [WH_1d()] function
#' @param newdata A list containing a vector indicating the new observation
#'   positions
#' @param ... Not used
#'
#' @returns An object of class \code{"WH_1d"} with additional components for model
#'   prediction.
#'
#' @examples
#' keep <- which(portfolios_mort[[1]]$ec > 0)
#' d <- portfolios_mort[[1]]$d[keep]
#' ec <- portfolios_mort[[1]]$ec[keep]
#'
#' WH_1d(d, ec) |> predict(newdata = 18:99) |> plot()
#'
#' @export
predict.WH_1d <- function(object, newdata = NULL, ...) {

  data <- as.numeric(names(object$y))
  full_data <- sort(union(data, newdata))
  ind <- order(c(data, setdiff(full_data, data)))

  n <- length(data)
  n_pred <- length(full_data)

  C <- diag(1L, n, n_pred)[, ind] # constraint matrix

  wt_pred <- c(t(C) %*% object$wt)
  W_pred <- diag(wt_pred) # extended weight matrix
  D_mat_pred <- build_D_mat(n_pred, object$q) # extended difference matrices
  P_pred <- object$lambda * crossprod(D_mat_pred) # extended penalization matrix

  Psi_pred <- (W_pred + P_pred) |> chol() |> chol2inv() # unconstrained variance / covariance matrix
  Psi_inv <- object$Psi |> chol() |> chol2inv()
  Psi_bis_inv <- (C %*% Psi_pred %*% t(C)) |> chol() |> chol2inv()

  A_pred <- Psi_pred %*% t(C) %*% Psi_bis_inv
  y_pred <- c(A_pred %*% object$y_hat)
  std_y_pred <- sqrt(colSums(t(A_pred) * (object$Psi %*% t(A_pred))))

  A_pred_2 <- Psi_pred %*% t(C) %*% Psi_inv
  y_pred_2 <- c(A_pred_2 %*% object$y_hat)
  std_y_pred_2 <- sqrt(colSums(t(A_pred_2) * (object$Psi %*% t(A_pred_2))))

  names(y_pred) <- names(std_y_pred) <-
    names(std_y_pred) <- names(std_y_pred_2) <- full_data

  object$y_pred <- y_pred
  object$std_y_pred <- std_y_pred
  object$y_pred_2 <- y_pred_2
  object$std_y_pred_2 <- std_y_pred_2

  return(object)
}

#' Predict new values using a fitted 2D WH model
#'
#' Extrapolate the model for new values of the covariates in a way that is
#' consistent with the fitted values
#'
#' @param object An object of class \code{"WH_2d"} returned by the [WH_2d()] function
#' @param newdata A list containing two vectors indicating the new observation
#'   positions
#' @param ... Not used
#'
#' @returns An object of class \code{"WH_2d"} with additional components for model
#'   prediction.
#'
#' @examples
#' keep_age <- which(rowSums(portfolios_LTC[[1]]$ec) > 1e2)
#' keep_duration <- which(colSums(portfolios_LTC[[1]]$ec) > 1e2)
#'
#' d  <- portfolios_LTC[[1]]$d[keep_age, keep_duration]
#' ec <- portfolios_LTC[[1]]$ec[keep_age, keep_duration]
#'
#' WH_2d(d, ec) |> predict(newdata = list(age = 50:99, duration = 0:19)) |> plot()
#'
#' @export
predict.WH_2d <- function(object, newdata = NULL, ...) {

  data <- dimnames(object$y) |> map(as.numeric)
  full_data <- map2(data, newdata, \(x,y) sort(union(x, y)))
  ind <- map2(data, full_data, \(x,y) order(c(x, setdiff(y, x))))

  n <- map_int(data, length)
  n_pred <- map_int(full_data, length)

  C <- map2(n, n_pred, \(x,y) diag(1L, x, y)) |>
    map2(ind, \(x,y) x[,y]) |>
    rev() |>
    reduce(kronecker) # constraint matrix

  wt_pred <- c(t(C) %*% c(object$wt))
  W_pred <- diag(wt_pred) # extended weight matrix
  D_mat_pred <- map2(n_pred, object$q, build_D_mat) # extended difference matrices
  P_pred <- object$lambda[[1]] * diag(n_pred[[2]]) %x% crossprod(D_mat_pred[[1]]) +
    object$lambda[[2]] * crossprod(D_mat_pred[[2]]) %x% diag(n_pred[[1]]) # extended penalization matrix

  Psi_pred <- (W_pred + P_pred) |> chol() |> chol2inv() # unconstrained variance / covariance matrix
  Psi_inv <- object$Psi |> chol() |> chol2inv()
  Psi_bis_inv <- (C %*% Psi_pred %*% t(C)) |> chol() |> chol2inv()

  A_pred <- Psi_pred %*% t(C) %*% Psi_bis_inv
  y_pred <- c(A_pred %*% c(object$y_hat))
  std_y_pred <- sqrt(colSums(t(A_pred) * (object$Psi %*% t(A_pred))))

  A_pred_2 <- Psi_pred %*% t(C) %*% Psi_inv
  y_pred_2 <- c(A_pred_2 %*% c(object$y_hat))
  std_y_pred_2 <- sqrt(colSums(t(A_pred_2) * (object$Psi %*% t(A_pred_2))))

  dim(y_pred) <- dim(std_y_pred) <-
    dim(y_pred_2) <- dim(std_y_pred_2) <- map_int(full_data, length) # set dimension for output matrices
  dimnames(y_pred) <- dimnames(std_y_pred) <-
    dimnames(y_pred_2) <- dimnames(std_y_pred_2) <- full_data # set names for output matrices

  object$y_pred <- y_pred
  object$std_y_pred <- std_y_pred
  object$y_pred_2 <- y_pred_2
  object$std_y_pred_2 <- std_y_pred_2

  return(object)
}

# Mise en forme----

#' Store WH model fit results in a data.frame
#'
#' @param object An object of class  \code{"WH_1d"} or \code{"WH_2d"} returned
#'   by one of the eponymous functions [WH_1d()] or [WH_2d()]
#' @param dim1 The (optional) name to be given to the first dimension
#' @param dim2 The (optional) name to be given to the second dimension
#'
#' @returns A data.frame gathering information about the fitted and predicted
#'   values, the model variance, residuals and effective degrees of freedom...
#'
#' @export
output_to_df <- function(object, dim1 = "x", dim2 = "t") {

  if ("y_pred" %in% names(object)) {
    object$y_hat <- object$y_pred
    object$std_y_hat = object$std_y_pred
  }

  df <- data.frame(y_hat = c(object$y_hat),
                   std_y_hat = c(object$std_y_hat))

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

    df[,c(dim1, intersect(
      c("ec", "d", "y", "y_hat", "std_y_hat", "wt", "res", "edf"),
      names(object)))]

  } else {

    x <- as.numeric(rownames(object$y))
    x_pred <- as.numeric(rownames(object$y_hat))

    t <- as.numeric(colnames(object$y))
    t_pred <- as.numeric(colnames(object$y_hat))

    df[[dim1]] <- rep(x_pred, times = length(t_pred))
    df[[dim2]] <- rep(t_pred, each = length(x_pred))

    df$y <- df$wt <- df$res <- df$edf <- NA_real_
    df$y[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$y)
    df$wt[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$wt)
    df$res[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$res)
    df$edf[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$edf)

    if("ec" %in% names(object) && "d" %in% names(object)) {

      df$ec <- df$d <- NA_real_

      df$ec[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$ec)
      df$d[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$d)
    }

    df[,c(dim1, dim2, intersect(
      c("ec", "d", "y", "y_hat", "std_y_hat", "wt", "res", "edf"),
      names(object)))]
  }
}

# Affichage----

#' Display of 1D WH object
#'
#' @param x An object of class \code{"WH_1d"} returned by the [WH_1d()] function
#' @param ... Not used
#'
#' @returns Invisibly returns \code{x}.
#'
#' @examples
#' keep <- which(portfolios_mort[[1]]$ec > 0)
#' d <- portfolios_mort[[1]]$d[keep]
#' ec <- portfolios_mort[[1]]$ec[keep]
#'
#' y <- log(d / ec)
#' y[d == 0] <- - 20
#' wt <- d
#'
#' WH_1d(d, ec) |> summary()
#' WH_1d(d, ec, method = "optim") |> summary()
#' WH_1d(d, ec, criterion = "GCV") |> summary()
#' WH_1d(y = y, wt = wt) |> summary()
#'
#' @export
print.WH_1d <- function(x, ...) {

  cat("\n")
  cat("An object fitted using the WH_1D function\n")
  cat(paste("Initial data contains", length(x$y), "data points\n"))
  cat(paste("Optimal smoothing parameter selected:", format(x$lambda, digits = 3), "\n"))
  cat(paste("Associated degrees of freedom:", format(sum(x$edf, digits = 3)), "\n"))
  invisible(x)
}

#' Display of 2D WH object
#'
#' @param x An object of class \code{"WH_2d"} returned by the [WH_2d()] function
#' @param ... Not used
#'
#' @returns Invisibly returns \code{x}.
#'
#' @examples
#' keep_age <- which(rowSums(portfolios_LTC[[1]]$ec) > 1e2)
#' keep_duration <- which(colSums(portfolios_LTC[[1]]$ec) > 1e2)
#'
#' d  <- portfolios_LTC[[1]]$d[keep_age, keep_duration]
#' ec <- portfolios_LTC[[1]]$ec[keep_age, keep_duration]
#'
#' WH_2d(d, ec) |> summary()
#'
#' @export
print.WH_2d <- function(x, ...) {

  cat("\n")
  cat("An object fitted using the WH_2D function\n")
  cat(paste("Initial data contains", prod(dim(x$y)), "data points\n"))
  cat(paste("Optimal smoothing parameters selected:", format(x$lambda, digits = 3), "\n"))
  cat(paste("Associated degrees of freedom:", format(sum(x$edf, digits = 3)), "\n"))
  invisible(x)
}

# Plots----

#' Plot 1D WH fit
#'
#' @param x An object of class \code{"WH_1d"} returned by the [WH_1d()] function
#' @param what What should be plotted. Should be one of \code{fit} (the
#'   default), \code{res} for residuals and \code{edf} for the effective degrees
#'   of freedom.
#' @param trans An (optional) transformation to be applied to the data. By
#'   default the identity function
#' @param ... Not used
#'
#' @returns A plot representing the desired element of the fit
#'
#' @examples
#' keep <- which(portfolios_mort[[1]]$ec > 0)
#' d <- portfolios_mort[[1]]$d[keep]
#' ec <- portfolios_mort[[1]]$ec[keep]
#'
#' WH_1d(d, ec) |> plot()
#' WH_1d(d, ec) |> plot("res")
#' WH_1d(d, ec) |> plot("edf")
#'
#' @export
plot.WH_1d <- function(x, what = "fit", trans, ...) {

  df <- output_to_df(x)
  if (missing(trans)) trans <- \(x) x

  switch(what,
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
           plot(df$x, df$res,
                xlab = "x", ylab = "deviance residuals", type = "b")
           graphics::abline(a = 0, b = 0, lty = 2, col = "blue")
         },
         edf = {
           plot(df$x, df$edf,
                xlab = "x", ylab = "degrees of freedom", type = "b")
           graphics::abline(a = 0, b = 0, lty = 2, col = "blue")
           graphics::abline(a = 1, b = 0, lty = 2, col = "blue")
         })
}

#' Plot 2D WH fit
#'
#' @param x An object of class \code{"WH_2d"} returned by the [WH_2d()] function
#' @param what What should be plotted (y_hat, std_y_hat, res, edf)
#' @param trans An (optional) transformation to be applied to the data
#' @param ... Not used
#'
#' @returns A plot representing the given element of the fit...
#'
#' @examples
#'
#' keep_age <- which(rowSums(portfolios_LTC[[1]]$ec) > 1e2)
#' keep_duration <- which(colSums(portfolios_LTC[[1]]$ec) > 1e2)
#'
#' d  <- portfolios_LTC[[1]]$d[keep_age, keep_duration]
#' ec <- portfolios_LTC[[1]]$ec[keep_age, keep_duration]
#'
#' WH_2d(d, ec) |> plot()
#' WH_2d(d, ec) |> plot("std_y_hat")
#'
#' @export
plot.WH_2d <- function(x, what = "y_hat", trans, ...) {

  df <- output_to_df(x)
  if (missing(trans)) trans <- \(x) x

  x <- unique(df$x)
  t <- unique(df$t)

  data <- matrix(c(df[[what]]), length(t), length(x), byrow = TRUE)

  graphics::contour(
    t, x, data,
    nlevels = 20,
    xlab = "x",
    ylab = "t",
    main = switch(
      what,
      y_hat = "log - hazard rate",
      std_y_hat = " standard deviation of log - hazard rate",
      edf = "degrees of freedom"))
}

# Regression----

## 1D----

#' Whittaker-Henderson Smoothing (Regression, fixed lambda)
#'
#' @param y Vector of observations
#' @param wt Optional vector of weights
#' @param lambda Smoothing parameter
#' @param p The number of eigenvectors to keep
#' @param q Order of penalization. Polynoms of degrees q - 1 are considered
#'   smooth and are therefore unpenalized
#'
#' @returns An object of class \code{"WH_1d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_1d_reg_fixed_lambda <- function(y, wt = rep(1, length(y)), lambda = 1e3, p = length(y), q = 2) {

  n <- length(y)
  SVD <- eigen_dec(n, q, p)

  X <- SVD$X
  Z <- SVD$Z
  U <- cbind(X, Z)

  s <- SVD$s
  s_tilde <- s[- seq_len(q)]

  tUWU <- t(U) %*% (wt * U)
  tUWy <- t(U) %*% (wt * y)

  Psi_chol <- tUWU
  diag(Psi_chol) <- diag(Psi_chol) + lambda * s
  Psi_chol <- Psi_chol |> chol()
  Psi <- Psi_chol |> chol2inv()

  gamma_hat <- c(Psi %*% tUWy)

  RESS <- sum(gamma_hat * s * gamma_hat)
  edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

  y_hat <- c(U %*% gamma_hat)
  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- wt * diag(Psi) # effective degrees of freedom by observation

  n_pos <- sum(wt != 0)
  dev <- sum(res * res) # residuals sum of squares
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- (p - q) * log(lambda) + sum(log(s_tilde))
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev + lambda * RESS - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(edf) <-
    names(wt) <- names(y) # set names for output vectors

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, p = p, q = q)
  class(out) <- "WH_1d"

  return(out)
}

#' Whittaker-Henderson Smoothing (Regression, optimize function)
#'
#' @inheritParams WH_1d_reg_fixed_lambda
#' @param criterion Criterion used to choose the smoothing parameter. One of
#'   "GCV" (default), "AIC" or "BIC".
#' @param lambda Initial smoothing parameter
#' @param verbose Should information about the optimization progress be
#'   displayed
#' @param accu_edf Tolerance for the convergence of the outer optimization
#'   procedure
#'
#' @returns An object of class \code{"WH_1d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_1d_reg_optim <- function(y, wt = rep(1, length(y)), p = length(y), q = 2,
                            criterion = "REML", lambda = 1e3,
                            verbose = FALSE, accu_edf = 1e-10) {

  n <- length(y)
  n_pos <- sum(wt != 0)
  SVD <- eigen_dec(n, q, p)

  X <- SVD$X
  Z <- SVD$Z
  U <- cbind(X, Z)

  s <- SVD$s
  s_tilde <- s[- seq_len(q)]

  tUWU <- t(U) %*% (wt * U)
  tUWy <- t(U) %*% (wt * y)

  WH_1d_reg_aux <- function(log_lambda) {

    lambda <- exp(log_lambda)
    if (verbose) cat("lambda : ", format(lambda, digits = 3), "\n")

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + lambda * s
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    gamma_hat <- c(Psi %*% tUWy)

    RESS <- sum(gamma_hat * s * gamma_hat)
    edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

    y_hat <- c(U %*% gamma_hat)

    res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
    dev <- sum(res * res) # residuals sum of squares
    sum_edf <- sum(edf_par) # effective degrees of freedom

    switch(criterion,
           AIC = dev + 2 * sum_edf,
           BIC = dev + log(n_pos) * sum_edf,
           GCV = n_pos * dev / (n_pos - sum_edf) ^ 2,
           REML = {
             tr_log_P <- (p - q) * log(lambda) + sum(log(s_tilde))
             tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

             dev + lambda * RESS - tr_log_P + tr_log_Psi
           })
  }

  lambda <- exp(stats::optimize(f = WH_1d_reg_aux, interval = 25 * c(- 1, 1), tol = accu_edf)$minimum)
  WH_1d_reg_fixed_lambda(y, wt, lambda, p, q)
}

#' Whittaker-Henderson Smoothing (Regression, Generalized Fellner-Schall update)
#'
#' @inheritParams WH_1d_reg_optim
#'
#' @returns An object of class \code{"WH_1d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_1d_reg_fs <- function(y, wt = rep(1, length(y)), p = length(y), q = 2,
                             lambda = 1e3, verbose = FALSE, accu_edf = 1e-10) {

  # Initialization
  n <- length(y)
  SVD <- eigen_dec(n, q, p)

  X <- SVD$X
  Z <- SVD$Z
  U <- cbind(X, Z)

  s <- SVD$s
  s_tilde <- s[- seq_len(q)]

  tUWU <- t(U) %*% (wt * U)
  tUWy <- t(U) %*% (wt * y)

  init <- TRUE

  # Loop
  while (init || cond_edf_random) {

    lambda <- if (init) lambda else (edf_random / RESS)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + lambda * s
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    gamma_hat <- c(Psi %*% tUWy)

    RESS <- sum(gamma_hat * s * gamma_hat)
    edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

    old_edf_random <- if (init) NA else edf_random
    edf_random <- sum(edf_par) - q
    if (verbose) cat("edf :", format(old_edf_random + q, digits = 3), "=>",
                     format(edf_random + q, digits = 3), "\n")
    cond_edf_random <- if (init) TRUE else {
      abs(edf_random - old_edf_random) > accu_edf * (old_edf_random + q)
    }
    init <- FALSE
  }

  y_hat <- c(U %*% gamma_hat)
  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- wt * diag(Psi) # effective degrees of freedom by observation

  n_pos <- sum(wt != 0)
  dev <- sum(res * res) # residuals sum of squares
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- (p - q) * log(lambda) + sum(log(s_tilde))
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev + lambda * RESS - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(edf) <-
    names(wt) <- names(y) # set names for output vectors

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, p = p, q = q)
  class(out) <- "WH_1d"

  return(out)
}

## 2D----

#' 2D Whittaker-Henderson Smoothing (Regression, fixed lambda)
#'
#' @param y Matrix of observations
#' @param wt Optional matrix of weights
#' @param lambda Vector of smoothing parameter
#' @param p The number of eigenvectors to keep on each dimension
#' @param q Matrix of orders of penalization. Polynoms of degrees q - 1 are considered
#'   smooth and are therefore unpenalized
#'
#' @returns An object of class \code{"WH_2d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_2d_reg_fixed_lambda <- function(y, wt = matrix(1, nrow = nrow(y), ncol = ncol(y)),
                                   lambda = c(1e3, 1e3), p = dim(y), q = c(2, 2)) {

  n <- dim(y)
  which_pos <- which(wt != 0)
  SVD <- pmap(list(n = n, q = q, p = p), eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()
  U <- cbind(X, Z)
  U_pos <- U[which_pos,]

  s_SVD <- map(SVD, "s")
  s_tilde_SVD <- map2(s_SVD, q, \(x, y) x[- seq_len(y)])
  s_tilde <- list(c(rep(s_tilde_SVD[[1]], q[[2]]),
                    rep(0, q[[1]] * (p[[2]] - q[[2]])),
                    rep(s_tilde_SVD[[1]], p[[2]] - q[[2]])),
                  c(rep(0, q[[2]] * (p[[1]] - q[[1]])),
                    rep(s_tilde_SVD[[2]], each = q[[1]]),
                    rep(s_tilde_SVD[[2]], each = p[[1]] - q[[1]])))
  s <- list(c(rep(0, prod(q)), s_tilde[[1]]),
            c(rep(0, prod(q)), s_tilde[[2]]))

  wt_pos <- c(wt)[which_pos]
  y_pos <- c(y)[which_pos]

  tUWU <- t(U_pos) %*% (wt_pos * U_pos)
  tUWy <- t(U_pos) %*% (wt_pos * y_pos)

  s_lambda <- map2(lambda, s, `*`)
  sum_s_lambda <- s_lambda |> do.call(what = `+`)

  Psi_chol <- tUWU
  diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
  Psi_chol <- Psi_chol |> chol()
  Psi <- Psi_chol |> chol2inv()

  gamma_hat <- c(Psi %*% tUWy)

  RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))
  edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter
  omega_j <- map(s_lambda, \(x) ifelse(x == 0, 0, x / sum_s_lambda))

  y_hat <- c(U %*% gamma_hat)
  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- c(wt) * diag(Psi) # effective degrees of freedom by observation / parameter

  n_pos <- sum(wt != 0)
  dev <- sum(res * res) # residuals sum of squares
  pen <- map2(lambda, RESS, `*`) |> do.call(what = `+`)
  dev_pen <- dev + pen
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- map2(lambda, s_tilde, `*`) |> do.call(what = `+`) |> log() |> sum()
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev + pen - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(edf) <-
    dim(wt) <- dim(y) # set dimensions for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(edf) <-
    dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, omega_j = omega_j, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, p = p, q = q)
  class(out) <- "WH_2d"

  return(out)
}

#' 2D Whittaker-Henderson Smoothing (Regression, optim function)
#'
#' @inheritParams WH_1d_reg_optim
#' @inheritParams WH_2d_reg_fixed_lambda
#' @param lambda Initial smoothing parameters
#'
#' @returns An object of class \code{"WH_2d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_2d_reg_optim <- function(y, wt = matrix(1, nrow = nrow(y), ncol = ncol(y)),
                            q = c(2, 2), p = dim(y), criterion = "REML",
                            lambda = c(1e3, 1e3), verbose = FALSE, accu_edf = 1e-10) {

  n <- dim(y)
  n_pos <- sum(wt != 0)
  which_pos <- which(wt != 0)
  SVD <- pmap(list(n = n, q = q, p = p), eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()
  U <- cbind(X, Z)
  U_pos <- U[which_pos,]

  s_SVD <- map(SVD, "s")
  s_tilde_SVD <- map2(s_SVD, q, \(x, y) x[- seq_len(y)])
  s_tilde <- list(c(rep(s_tilde_SVD[[1]], q[[2]]),
                    rep(0, q[[1]] * (p[[2]] - q[[2]])),
                    rep(s_tilde_SVD[[1]], p[[2]] - q[[2]])),
                  c(rep(0, q[[2]] * (p[[1]] - q[[1]])),
                    rep(s_tilde_SVD[[2]], each = q[[1]]),
                    rep(s_tilde_SVD[[2]], each = p[[1]] - q[[1]])))
  s <- list(c(rep(0, prod(q)), s_tilde[[1]]),
            c(rep(0, prod(q)), s_tilde[[2]]))

  wt_pos <- c(wt)[which_pos]
  y_pos <- c(y)[which_pos]

  tUWU <- t(U_pos) %*% (wt_pos * U_pos)
  tUWy <- t(U_pos) %*% (wt_pos * y_pos)

  WH_2d_reg_aux <- function(log_lambda) {

    lambda <- exp(log_lambda)
    if (verbose) cat("lambda : ", format(lambda, digits = 3), "\n")

    s_lambda <- map2(lambda, s, `*`)
    sum_s_lambda <- s_lambda |> do.call(what = `+`)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    gamma_hat <- c(Psi %*% tUWy)

    RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))
    edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

    y_hat <- c(U %*% gamma_hat)

    res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
    dev <- sum(res * res) # residuals sum of squares
    sum_edf <- sum(edf_par) # effective degrees of freedom

    switch(criterion,
           AIC = dev + 2 * sum_edf,
           BIC = dev + log(prod(n_pos)) * sum_edf,
           GCV = prod(n_pos) * dev / (prod(n_pos) - sum_edf) ^ 2,
           REML = {
             tr_log_P <- map2(lambda, s_tilde, `*`) |> do.call(what = `+`) |> log() |> sum()
             tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

             dev + lambda[[1]] * RESS[[1]] + lambda[[2]] * RESS[[2]] - tr_log_P + tr_log_Psi
           })
  }

  lambda <- exp(stats::optim(par = log(lambda), fn = WH_2d_reg_aux, control = list(reltol = accu_edf))$par)
  WH_2d_reg_fixed_lambda(y, wt, lambda, p, q)
}

#' 2D Whittaker-Henderson Smoothing (Regression, Fellner-Schall update)
#'
#' @inheritParams WH_2d_reg_optim
#'
#' @returns An object of class \code{"WH_2d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_2d_reg_fs <- function(y, wt = matrix(1, nrow = nrow(y), ncol = ncol(y)),
                         q = c(2, 2), p = dim(y), lambda = c(1e3, 1e3),
                         verbose = FALSE, accu_edf = 1e-10) {

  # Initialization
  n <- dim(y)
  which_pos <- which(wt != 0)
  SVD <- pmap(list(n = n, q = q, p = p), eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()
  U <- cbind(X, Z)
  U_pos <- U[which_pos,]

  s_SVD <- map(SVD, "s")
  s_tilde_SVD <- map2(s_SVD, q, \(x, y) x[- seq_len(y)])
  s_tilde <- list(c(rep(s_tilde_SVD[[1]], q[[2]]),
                    rep(0, q[[1]] * (p[[2]] - q[[2]])),
                    rep(s_tilde_SVD[[1]], p[[2]] - q[[2]])),
                  c(rep(0, q[[2]] * (p[[1]] - q[[1]])),
                    rep(s_tilde_SVD[[2]], each = q[[1]]),
                    rep(s_tilde_SVD[[2]], each = p[[1]] - q[[1]])))
  s <- list(c(rep(0, prod(q)), s_tilde[[1]]),
            c(rep(0, prod(q)), s_tilde[[2]]))

  wt_pos <- c(wt)[which_pos]
  y_pos <- c(y)[which_pos]

  tUWU <- t(U_pos) %*% (wt_pos * U_pos)
  tUWy <- t(U_pos) %*% (wt_pos * y_pos)

  init <- TRUE

  # Loop
  while (init || cond_edf_random) {

    if (!init) lambda <- edf_random / RESS

    s_lambda <- map2(lambda, s, `*`)
    sum_s_lambda <- s_lambda |> do.call(what = `+`)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    gamma_hat <- c(Psi %*% tUWy)

    RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))
    edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter
    omega_j <- map(s_lambda, \(x) ifelse(x == 0, 0, x / sum_s_lambda))

    old_edf_random <- if (init) NA else edf_random
    edf_random <- map_dbl(omega_j, \(x) sum(x * edf_par))
    if (verbose) cat("edf :", format(old_edf_random + q, digits = 3),
                     "=>", format(edf_random + q, digits = 3), "\n")
    cond_edf_random <- if (init) TRUE else{
      any(abs(edf_random - old_edf_random) > accu_edf * (old_edf_random + q))
    }
    init <- FALSE
  }

  y_hat <- c(U %*% gamma_hat)
  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- sqrt(wt) * (y - y_hat) # (weighted) residuals
  edf <- c(wt) * diag(Psi) # effective degrees of freedom by observation / parameter

  n_pos <- sum(wt != 0)
  dev <- sum(res * res) # residuals sum of squares
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- map2(lambda, s_tilde, `*`) |> do.call(what = `+`) |> log() |> sum()
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev + lambda[[1]] * RESS[[1]] + lambda[[2]] * RESS[[2]] - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(edf) <-
    dim(wt) <- dim(y) # set dimensions for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(edf) <-
    dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, omega_j = omega_j, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, p = p, q = q)
  class(out) <- "WH_2d"

  return(out)
}

# Maximum Likelihood----

## 1D----

#' Whittaker-Henderson Smoothing (Maximum Likelihood, fixed lambda)
#'
#' @inheritParams WH_1d_reg_fixed_lambda
#' @inheritParams WH_1d_reg_optim
#' @param d Vector of observed events
#' @param ec Vector of central exposure
#' @param accu_dev Tolerance for the convergence of the optimization procedure
#'
#' @returns An object of class \code{"WH_1d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_1d_ml_fixed_lambda <- function(d, ec, lambda = 1e3, p = length(d), q = 2,
                                  verbose = FALSE, accu_dev = 1e-12) {

  # Initialization
  n <- length(d)
  sum_d <- sum(d)
  off <- log(pmax(ec, 1e-4))
  y <- ifelse(d == 0, NA, log(d)) - off

  SVD <- eigen_dec(n, q, p)

  X <- SVD$X
  Z <- SVD$Z
  U <- cbind(X, Z)

  s <- SVD$s
  s_tilde <- s[- seq_len(q)]

  y_hat <- log(pmax(d, 1e-8)) - off
  new_wt <- exp(y_hat + off)
  # gamma_hat <- c(t(U) %*% y_hat)
  # RESS <- sum(gamma_hat * s * gamma_hat)
  # dev_pen <- compute_deviance(d, new_wt) + lambda * RESS
  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    # update of parameters, working vector and weight matrix
    wt <- new_wt
    z <- y_hat + d / wt - 1

    tUWU <- t(U) %*% (wt * U)
    tUWz <- t(U) %*% (wt * z)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + lambda * s
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    gamma_hat <- c(Psi %*% tUWz) # fitted value

    y_hat <- c(U %*% gamma_hat)
    new_wt <- exp(y_hat + off)

    # update of convergence check
    old_dev_pen <- dev_pen
    RESS <- sum(gamma_hat * s * gamma_hat)
    dev <- compute_deviance(d, new_wt)
    dev_pen <- dev + lambda * RESS
    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- (old_dev_pen - dev_pen) > accu_dev * sum_d
  }

  edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- compute_res_deviance(d, new_wt) # (weighted) residuals
  edf <- wt * diag(Psi) # effective degrees of freedom by observation / parameter

  n_pos <- sum(wt != 0)
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- (p - q) * log(lambda) + sum(log(s_tilde))
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev_pen - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(edf) <-
    names(wt) <- names(z) <- names(y) # set names for output vectors

  out <- list(d = d, ec = ec, y = y, wt = wt, z = z, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, q = q)
  class(out) <- "WH_1d"

  return(out)
}

#' Whittaker-Henderson Smoothing (Maximum Likelihood, optimize function)
#'
#' @inheritParams WH_1d_reg_optim
#' @inheritParams WH_1d_ml_fixed_lambda
#'
#' @returns An object of class \code{"WH_1d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_1d_ml_optim <- function(d, ec, p = length(d), q = 2, criterion = "REML", lambda = 1e3,
                             verbose = FALSE, accu_edf = 1e-10, accu_dev = 1e-12) {

  # Initialization
  n <- length(d)
  sum_d <- sum(d)
  off <- log(pmax(ec, 1e-4))
  y <- ifelse(d == 0, NA, log(d)) - off

  SVD <- eigen_dec(n, q, p)

  X <- SVD$X
  Z <- SVD$Z
  U <- cbind(X, Z)

  s <- SVD$s
  s_tilde <- s[- seq_len(q)]

  WH_1d_ml_aux <- function(log_lambda) {

    lambda <- exp(log_lambda)
    if (verbose) cat("lambda : ", format(lambda, digits = 3), "\n")

    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
    # gamma_hat <- c(t(U) %*% y_hat)
    # RESS <- sum(gamma_hat * s * gamma_hat)
    # dev_pen <- compute_deviance(d, new_wt) + lambda * RESS
    dev_pen <- Inf
    cond_dev_pen <- TRUE

    # Loop
    while (cond_dev_pen) {

      # update of working vector and weight matrix
      wt <- new_wt
      z <- y_hat + d / wt - 1

      tUWU <- t(U) %*% (wt * U)
      tUWz <- t(U) %*% (wt * z)

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + lambda * s
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      gamma_hat <- c(Psi %*% tUWz) # fitted value

      RESS <- sum(gamma_hat * s * gamma_hat)

      y_hat <- c(U %*% gamma_hat)
      new_wt <- exp(y_hat + off)

      # update of convergence check
      old_dev_pen <- dev_pen
      dev <- compute_deviance(d, new_wt)
      dev_pen <- dev + lambda * RESS
      if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                       "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
      cond_dev_pen <- (old_dev_pen - dev_pen) > accu_dev * sum_d
    }

    n_pos <- sum(wt != 0)
    sum_edf <- sum(t(Psi) * tUWU) # effective degrees of freedom

    switch(criterion,
           AIC = dev + 2 * sum_edf,
           BIC = dev + log(n_pos) * sum_edf,
           GCV = n_pos * dev / (n_pos - sum_edf) ^ 2,
           REML = {
             tr_log_P <- (p - q) * log(lambda) + sum(log(s_tilde))
             tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

             dev_pen - tr_log_P + tr_log_Psi
           })
  }

  lambda <- exp(stats::optimize(f = WH_1d_ml_aux, interval = 25 * c(- 1, 1), tol = accu_edf)$minimum)
  WH_1d_ml_fixed_lambda(d, ec, lambda, p, q)
}

#' Whittaker-Henderson Smoothing (Maximum Likelihood, Generalized Fellner-Schall update)
#'
#' @inheritParams WH_1d_ml_optim
#'
#' @returns An object of class \code{"WH_1d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_1d_ml_fs <- function(d, ec, p = length(d), q = 2,
                        lambda = 1e3, verbose = FALSE, accu_edf = 1e-10, accu_dev = 1e-12) {

  # Initialization
  n <- length(d)
  sum_d <- sum(d)
  off <- log(pmax(ec, 1e-4))
  y <- ifelse(d == 0, NA, log(d)) - off

  SVD <- eigen_dec(n, q, p)

  X <- SVD$X
  Z <- SVD$Z
  U <- cbind(X, Z)

  s <- SVD$s
  s_tilde <- s[- seq_len(q)]

  y_hat <- log(pmax(d, 1e-8)) - off
  new_wt <- exp(y_hat + off)

  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    # update of working vector and weight matrix
    wt <- new_wt
    z <- y_hat + d / wt - 1

    tUWU <- t(U) %*% (wt * U)
    tUWz <- t(U) %*% (wt * z)

    init_lambda <- TRUE

    # Loop
    while (init_lambda || cond_edf_random) {

      if (!init_lambda) lambda <- edf_random / RESS

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + lambda * s
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      gamma_hat <- c(Psi %*% tUWz) # fitted value

      RESS <- sum(gamma_hat * s * gamma_hat)
      edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

      old_edf_random <- if (init_lambda) NA else edf_random
      edf_random <- sum(edf_par) - q
      if (verbose) cat("edf :", format(old_edf_random + q, digits = 3),
                       "=>", format(edf_random + q, digits = 3), "\n")
      cond_edf_random <- if (init_lambda) TRUE else {
        abs(edf_random - old_edf_random) > accu_edf * (old_edf_random + q)
      }
      init_lambda <- FALSE
    }

    y_hat <- c(U %*% gamma_hat)
    new_wt <- exp(y_hat + off)

    # update of convergence check
    old_dev_pen <- dev_pen
    dev <- compute_deviance(d, new_wt)
    dev_pen <- dev + lambda * RESS
    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- (old_dev_pen - dev_pen) > accu_dev * sum_d
  }

  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- compute_res_deviance(d, new_wt) # (weighted) residuals
  edf <- wt * diag(Psi) # effective degrees of freedom by observation / parameter

  n_pos <- sum(wt != 0)
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- (p - q) * log(lambda) + sum(log(s_tilde))
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev_pen - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(edf) <-
    names(wt) <- names(z) <- names(y) # set names for output vectors

  out <- list(d = d, ec = ec, y = y, wt = wt, z = z, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, q = q)
  class(out) <- "WH_1d"

  return(out)
}

## 2D----

#' 2D Whittaker-Henderson Smoothing (Maximum Likelihood, fixed lambda)
#'
#' @inheritParams WH_2d_reg_fixed_lambda
#' @inheritParams WH_1d_ml_fixed_lambda
#'
#' @returns An object of class \code{"WH_2d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_2d_ml_fixed_lambda <- function(d, ec, lambda = c(1e3, 1e3), p = dim(d), q = c(2, 2),
                                  verbose = FALSE, accu_dev = 1e-12) {

  # Initialization
  n <- dim(d)
  which_pos <- which(ec != 0)
  sum_d <- sum(d)
  off <- log(pmax(ec, 1e-4))
  y <- ifelse(d == 0, NA, log(d)) - off

  SVD <- pmap(list(n = n, q = q, p = p), eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()
  U <- cbind(X, Z)
  U_pos <- U[which_pos,]

  s_SVD <- map(SVD, "s")
  s_tilde_SVD <- map2(s_SVD, q, \(x, y) x[- seq_len(y)])
  s_tilde <- list(c(rep(s_tilde_SVD[[1]], q[[2]]),
                    rep(0, q[[1]] * (p[[2]] - q[[2]])),
                    rep(s_tilde_SVD[[1]], p[[2]] - q[[2]])),
                  c(rep(0, q[[2]] * (p[[1]] - q[[1]])),
                    rep(s_tilde_SVD[[2]], each = q[[1]]),
                    rep(s_tilde_SVD[[2]], each = p[[1]] - q[[1]])))
  s <- list(c(rep(0, prod(q)), s_tilde[[1]]),
            c(rep(0, prod(q)), s_tilde[[2]]))

  s_lambda <- map2(lambda, s, `*`)
  sum_s_lambda <- s_lambda |> do.call(what = `+`)

  y_hat <- log(pmax(d, 1e-8)) - off
  new_wt <- exp(y_hat + off)
  # gamma_hat <- c(t(U) %*% c(y_hat))
  # RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))
  # dev_pen <- compute_deviance(d, new_wt) + map2(lambda, RESS, `*`) |> do.call(what = `+`)
  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    # update of parameters, working vector and weight matrix
    wt <- new_wt
    z <- y_hat + d / wt - 1

    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]

    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    gamma_hat <- c(Psi %*% tUWz) # fitted value

    RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))

    y_hat <- c(U %*% gamma_hat)
    new_wt <- exp(y_hat + off)

    # update of convergence check
    old_dev_pen <- dev_pen
    dev <- compute_deviance(d, new_wt)
    dev_pen <- dev + map2(lambda, RESS, `*`) |> do.call(what = `+`)
    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- (old_dev_pen - dev_pen) > accu_dev * sum_d
  }

  edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter

  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- compute_res_deviance(d, new_wt) # (weighted) residuals
  edf <- wt * diag(Psi) # effective degrees of freedom by observation / parameter

  n_pos <- sum(wt != 0)
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- map2(lambda, s_tilde, `*`) |> do.call(what = `+`) |> log() |> sum()
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev_pen - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(edf) <-
    dim(wt) <- dim(y) # set dimensions for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(edf) <-
    dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(d = d, ec = ec, y = y, wt = wt, z = z, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, q = q)
  class(out) <- "WH_2d"

  return(out)
}

#' 2D Whittaker-Henderson Smoothing (Maximum Likelihood, optim function)
#'
#' @inheritParams WH_2d_ml_fixed_lambda
#' @inheritParams WH_2d_reg_optim
#'
#' @returns An object of class \code{"WH_2d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_2d_ml_optim <- function(d, ec, q = c(2, 2), p = dim(d), criterion = "REML", lambda = c(1e3, 1e3),
                             verbose = FALSE, accu_edf = 1e-10, accu_dev = 1e-12) {

  n <- dim(d)
  which_pos <- which(ec != 0)
  sum_d <- sum(d)
  off <- log(pmax(ec, 1e-4))
  y <- ifelse(d == 0, NA, log(d)) - off

  SVD <- pmap(list(n = n, q = q, p = p), eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()
  U <- cbind(X, Z)
  U_pos <- U[which_pos,]

  s_SVD <- map(SVD, "s")
  s_tilde_SVD <- map2(s_SVD, q, \(x, y) x[- seq_len(y)])
  s_tilde <- list(c(rep(s_tilde_SVD[[1]], q[[2]]),
                    rep(0, q[[1]] * (p[[2]] - q[[2]])),
                    rep(s_tilde_SVD[[1]], p[[2]] - q[[2]])),
                  c(rep(0, q[[2]] * (p[[1]] - q[[1]])),
                    rep(s_tilde_SVD[[2]], each = q[[1]]),
                    rep(s_tilde_SVD[[2]], each = p[[1]] - q[[1]])))
  s <- list(c(rep(0, prod(q)), s_tilde[[1]]),
            c(rep(0, prod(q)), s_tilde[[2]]))

  WH_2d_ml_aux <- function(log_lambda) {

    lambda <- exp(log_lambda)
    if (verbose) cat("lambda : ", format(lambda, digits = 3), "\n")

    s_lambda <- map2(lambda, s, `*`)
    sum_s_lambda <- s_lambda |> do.call(what = `+`)

    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
    # gamma_hat <- c(t(U) %*% c(y_hat))
    # RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))
    # dev_pen <- compute_deviance(d, new_wt) + map2(lambda, RESS, `*`) |> do.call(what = `+`)
    dev_pen <- Inf
    cond_dev_pen <- TRUE

    # Loop
    while (cond_dev_pen) {

      # update of parameters, working vector and weight matrix
      wt <- new_wt
      z <- y_hat + d / wt - 1

      wt_pos <- c(wt)[which_pos]
      z_pos <- c(z)[which_pos]

      tUWU <- t(U_pos) %*% (wt_pos * U_pos)
      tUWz <- t(U_pos) %*% (wt_pos * z_pos)

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      gamma_hat <- c(Psi %*% tUWz) # fitted value

      RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))

      y_hat <- c(U %*% gamma_hat)
      new_wt <- exp(y_hat + off)

      # update of convergence check
      old_dev_pen <- dev_pen
      dev <- compute_deviance(d, new_wt)
      dev_pen <- dev + map2(lambda, RESS, `*`) |> do.call(what = `+`)
      if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                       "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
      cond_dev_pen <- (old_dev_pen - dev_pen) > accu_dev * sum_d
    }

    n_pos <- sum(wt != 0)
    sum_edf <- sum(t(Psi) * tUWU) # effective degrees of freedom

    switch(criterion,
           AIC = dev + 2 * sum_edf,
           BIC = dev + log(prod(n_pos)) * sum_edf,
           GCV = prod(n_pos) * dev / (prod(n_pos) - sum_edf) ^ 2,
           REML = {
             tr_log_P <- map2(lambda, s_tilde, `*`) |> do.call(what = `+`) |> log() |> sum()
             tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

             dev_pen - tr_log_P + tr_log_Psi
           })
  }

  lambda <- exp(stats::optim(par = log(lambda),
                      fn = WH_2d_ml_aux,
                      control = list(reltol = accu_edf))$par)
  WH_2d_ml_fixed_lambda(d, ec, lambda, p, q)
}

#' 2D Whittaker-Henderson Smoothing (Maximum Likelihood, Generalized Fellner-Schall update)
#'
#' @inheritParams WH_2d_ml_optim
#'
#' @returns An object of class \code{"WH_2d"} i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
WH_2d_ml_fs <- function(d, ec, q = c(2, 2), p = dim(d), verbose = FALSE,
                          lambda = c(1e3, 1e3), accu_edf = 1e-10, accu_dev = 1e-12) {

  # Initialization
  n <- dim(d)
  which_pos <- which(ec != 0)
  sum_d <- sum(d)
  off <- log(pmax(ec, 1e-4))
  y <- ifelse(d == 0, NA, log(d)) - off

  SVD <- pmap(list(n = n, q = q, p = p), eigen_dec)

  X_SVD <- map(SVD, "X")
  Z_SVD <- map(SVD, "Z")
  X <- list(XX = list(X_SVD[[1]], X_SVD[[2]])) |>
    compute_XZ_mat()
  Z <- list(ZX = list(Z_SVD[[1]], X_SVD[[2]]),
            XZ = list(X_SVD[[1]], Z_SVD[[2]]),
            ZZ = list(Z_SVD[[1]], Z_SVD[[2]])) |>
    compute_XZ_mat()
  U <- cbind(X, Z)
  U_pos <- U[which_pos,]

  s_SVD <- map(SVD, "s")
  s_tilde_SVD <- map2(s_SVD, q, \(x, y) x[- seq_len(y)])
  s_tilde <- list(c(rep(s_tilde_SVD[[1]], q[[2]]),
                    rep(0, q[[1]] * (p[[2]] - q[[2]])),
                    rep(s_tilde_SVD[[1]], p[[2]] - q[[2]])),
                  c(rep(0, q[[2]] * (p[[1]] - q[[1]])),
                    rep(s_tilde_SVD[[2]], each = q[[1]]),
                    rep(s_tilde_SVD[[2]], each = p[[1]] - q[[1]])))
  s <- list(c(rep(0, prod(q)), s_tilde[[1]]),
            c(rep(0, prod(q)), s_tilde[[2]]))

  y_hat <- log(pmax(d, 1e-8)) - off
  new_wt <- exp(y_hat + off)
  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    # update of parameters, working vector and weight matrix
    wt <- new_wt
    z <- y_hat + d / wt - 1

    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]

    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

    init_lambda <- TRUE

    # Loop
    while (init_lambda || cond_edf_random) {

      if (!init_lambda) lambda <- edf_random / RESS

      s_lambda <- map2(lambda, s, `*`)
      sum_s_lambda <- s_lambda |> do.call(what = `+`)

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      gamma_hat <- c(Psi %*% tUWz) # fitted value

      RESS <- map_dbl(s, \(x) sum(gamma_hat * x * gamma_hat))
      edf_par <- colSums(t(Psi) * tUWU) # effective degrees of freedom by parameter
      omega_j <- map(s_lambda, \(x) ifelse(x == 0, 0, x / sum_s_lambda))

      old_edf_random <- if (init_lambda) NA else edf_random
      edf_random <- map_dbl(omega_j, \(x) sum(x * edf_par))
      if (verbose) cat("edf :", format(old_edf_random + q, digits = 3),
                       "=>", format(edf_random + q, digits = 3), "\n")
      cond_edf_random <- if (init_lambda) TRUE else any(abs(edf_random - old_edf_random) > accu_edf * (old_edf_random + q))
      init_lambda <- FALSE
    }

    y_hat <- c(U %*% gamma_hat)
    new_wt <- exp(y_hat + off)

    # update of convergence check
    old_dev_pen <- dev_pen
    dev <- compute_deviance(d, new_wt)
    dev_pen <- dev + map2(lambda, RESS, `*`) |> do.call(what = `+`)
    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- (old_dev_pen - dev_pen) > accu_dev * sum_d
  }

  Psi <- U %*% Psi %*% t(U)
  std_y_hat <- sqrt(diag(Psi)) # standard deviation of fit

  res <- compute_res_deviance(d, new_wt) # (weighted) residuals
  edf <- wt * diag(Psi) # effective degrees of freedom by observation / parameter

  n_pos <- sum(wt != 0)
  sum_edf <- sum(edf) # effective degrees of freedom

  tr_log_P <- map2(lambda, s_tilde, `*`) |> do.call(what = `+`) |> log() |> sum()
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  AIC <- dev + 2 * sum_edf
  BIC <- dev + log(n_pos) * sum_edf
  GCV <- n_pos * dev / (n_pos - sum_edf) ^ 2
  REML <- dev_pen - tr_log_P + tr_log_Psi

  diagnosis <- data.frame(sum_edf = sum_edf, AIC = AIC, BIC = BIC, GCV = GCV, REML = REML)

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(edf) <-
    dim(wt) <- dim(y) # set dimensions for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(edf) <-
    dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(d = d, ec = ec, y = y, wt = wt, z = z, y_hat = y_hat, std_y_hat = std_y_hat,
              res = res, edf = edf, edf_par = edf_par, diagnosis = diagnosis,
              Psi = Psi, lambda = lambda, q = q)
  class(out) <- "WH_2d"

  return(out)
}


