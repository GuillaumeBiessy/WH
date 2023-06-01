# Wrappers functions----

#' 1D Whittaker-Henderson Smoothing
#'
#' Main package function to apply Whittaker-Henderson smoothing in a
#' unidimensional survival analysis framework. It takes as input a vector of
#' observed events and a vector of associated central exposure, both
#' depending on a single covariate, and build a smooth version of the log-hazard
#' rate. Smoothing parameters may be supplied or automatically chosen according
#' to an adequate criterion such as `"REML"` (the default), `"AIC"`, `"BIC"` or
#' `"GCV"`. Whittaker-Henderson may be applied in a full maximum likelihood
#' framework (the default) or an approximate gaussian framework.
#'
#' @param d Vector of observed events, with named elements.
#' @param ec Vector of central exposure. The central exposure corresponds to the
#'   sum of the exposure period over the insured population. An individual
#'   experiencing an event of interest during the year will no longer be exposed
#'   afterward and the exposure should be reduced accordingly.
#' @param lambda Smoothing parameter. If missing, an optimization procedure will
#'   be used to find the optimal smoothing parameter. If supplied, no optimal
#'   smoothing parameter search will take place unless the `criterion` argument
#'   is also supplied, in which case `lambda` will be used as the starting
#'   parameter for the optimization procedure.
#' @param criterion Criterion to be used for the selection of the optimal
#'   smoothing parameter. Default is `"REML"` which stands for restricted
#'   maximum likelihood. Other options include `"AIC"`, `"BIC"` and `"GCV"`.
#' @param method Method to be used to find the optimal smoothing parameter.
#'   Default to `"fixed_lambda"` if `lambda` is supplied, meaning no
#'   optimization is performed. Otherwise, default to `"perf"` which means the
#'   performance iteration method is used. The alternative `"outer"` method is
#'   guaranteed to converge but slower. Both those methods rely on the
#'   `optimize` function from package `stats`.
#' @param q Order of penalization. Polynoms of degrees `q - 1` are considered
#'   smooth and are therefore unpenalized. Should be left to the default of `2`
#'   for most practical applications.
#' @param framework Default framework is `"ml"` which stands for maximum
#'   likelihood unless the `y` argument is also provided, in which case an
#'   `"reg"` or regression framework is used. The regression framework is an
#'   approximate gaussian framework is used instead.
#' @param y Optional vector of observations whose elements should be named. Used
#'   only in the regression framework and even in this case will be
#'   automatically computed from the `d` and `ec` arguments if those are
#'   supplied. May be useful when using Whittaker-Henderson smoothing outside of
#'   the survival analysis framework.
#' @param wt Optional vector of weights. As for the observation vector `y`, used
#'   only in the regression framework and even in this case will be
#'   automatically computed if the `d` argument is supplied. May be useful when
#'   using Whittaker-Henderson smoothing outside of the survival analysis
#'   framework.
#' @param quiet Should messages and warnings be displayed ? Default to `FALSE`,
#'   may be set to `TRUE` if numerous calls to this function are made.
#' @param ... Additional parameters passed to the smoothing function called.
#'
#' @returns An object of class `WH_1d` i.e. a list containing :
#' * `d` The inputed vector of observed events (if supplied as input)
#' * `ec` The inputed vector of central exposure (if supplied as input)
#' * `y` The observation vector (either supplied or computed as y = log(d) - log(ec))
#' * `wt` The inputed vector of weights (either supplied or computed as `d`)
#' * `y_hat` The vector of fitted value
#' * `std_y_hat` The vector of standard deviation associated with the fitted value
#' * `res` The vector of model deviance residuals
#' * `edf` The vector of effective degrees of freedom associated with each observation
#' * `edf_par` The vector of effective degrees of freedom associated with each eigenvector
#' * `diagnosis` A data.frame with one line containing the effective degrees of freedom of the model, the deviance of the fit as well as the AIC, BIC, GCV and REML criteria
#' * `Psi` The variance-covariance matrix associated with the fit, required for
#'   extrapolation.
#' * `lambda` The smoothing parameter used, either supplied or computed.
#' * `p` The number of components kept after rank reduction is used (it should not be in the one-dimension case).
#' * `q` The supplied order for the penalization.
#'
#' @examples
#' keep <- which(portfolio_mort$ec > 0)
#' d <- portfolio_mort$d[keep]
#' ec <- portfolio_mort$ec[keep]
#'
#' y <- log(d / ec)
#' y[d == 0] <- - 20
#' wt <- d
#'
#' # Maximum likelihood
#' WH_1d(d, ec, lambda = 1e2)
#' WH_1d(d, ec) # default performance iteration method based on the optimize function
#'
#' testthat::expect_equal(WH_1d(d, ec),
#'                        WH_1d(d, ec, method = "outer"),
#'                        tolerance = 1e-1)
#' # performance iteration method is approximate in maximum likelihood framework
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
#' testthat::expect_equal(WH_1d(y = y, wt = wt),
#'                        WH_1d(y = y, wt = wt, method = "outer"),
#'                        tolerance = 1e-5)
#' # performance iteration method is exact in regression framework
#'
#' @export
WH_1d <- function(d, ec, lambda, criterion, method, q = 2, framework, y, wt, quiet = FALSE, ...) {

  if (missing(framework)) framework <- if (!missing(y)) "reg" else "ml"
  if (length(framework) != 1 || !(framework %in% c("reg", "ml"))) {
    stop("framework should be one of ml (the default) or reg")
  }
  if (framework == "reg") {
    if (missing(y)) {
      if (!missing(d) && !missing(ec)) {
        if (!is.numeric(d)) stop("d should be a numeric vector")
        if (!is.numeric(ec)) stop("ec should be a numeric vector")
        if (length(d) != length(ec)) stop("Length of d and ec must match")
        if (!quiet) message("Computing y as log(d / ec)")
        y <- log(d / ec)
        y[d == 0] <- - 20
      } else {
        stop("Either y or both d and ec required in the regression framework")
      }
    } else {
      if (!is.numeric(y)) stop("y should be a numeric vector")
    }
    if (missing(wt)) {
      if (!missing(d)) {
        if (!quiet) message("Using d as weights")
        wt <- d
      } else {
        stop("Either wt or d needs to be supplied in the regression framework")
      }
    } else {
      if (!is.numeric(wt)) stop("wt should be a numeric vector")
    }
    if (length(y) != length(wt)) stop("Length of y and wt must match")
    if (is.null(names(y))) {
      y_names <- seq_along(y)
      if (!quiet) warning("No names found for y, setting names to: ", paste(range(y_names), collapse = " - "))
    }
  } else {
    if (missing(d)) stop("d argument required in the maximum likelihood framework")
    if (missing(ec)) stop("ec argument required in the maximum likelihood framework")
    if (!is.numeric(d)) stop("d should be a numeric vector")
    if (!is.numeric(ec)) stop("ec should be a numeric vector")
    if (length(d) != length(ec)) stop("Length of d and ec must match")
    if (is.null(names(d))) {
      d_names <- seq_along(d)
      if (!quiet) warning("No names found for d, setting names to: ", paste(range(d_names), collapse = " - "))
    }
  }

  if (missing(criterion)) criterion <- "REML"
  criteria <- c("REML", "AIC", "BIC", "GCV")
  if (length(criterion) != 1 || !(criterion %in% criteria)) stop(
    "criterion should be one of ", paste(criteria, collapse = ", "))

  if (missing(method)) {
    if (!missing(lambda)) {
      if (length(lambda) != 1) stop("smoothing parameter lambda should be of length 1")
      if (!is.numeric(lambda) || lambda <= 0) stop(
        "smoothing parameter lambda should be a strictly positive number")
      method <- "fixed_lambda"
    } else {
      if (!quiet) message("Using outer iteration / Brent method")
      method <- "outer"
    }
  }
  methods <- c("fixed_lambda", "perf", "outer")
  if (length("method") != 1 || !(method %in% methods)) stop(
    "method should be one of ", paste(methods, collapse = ", "))

  if (missing(lambda)) {
    lambda <- 1e3
  } else {
    if (method != "fixed_lambda" && !quiet) message("Both method and lambda specified, lambda used as starting parameter")
  }

  if (!is.numeric(q) || length(q) != 1 || q <= 0 ||
      (abs(q - round(q)) > .Machine$double.eps ^ 0.5)) stop(
        "q should be a strictly positive integer"
      )
  if (q > 3 && !quiet) warning("Differences of order q > 3 have no meaningful interpretation.\n",
                               "Consider using a lower value for q")

  reg <- (framework == "reg")
  what <- paste("WH_1d", method, sep = "_")
  args <- (if (reg) list(y = y, wt = wt) else list(d = d, ec = ec)) |>
    c(list(q = q, lambda = lambda, reg = reg),
      if (method %in% c("perf", "outer")) list(criterion = criterion),
      list(...))

  out <- do.call(what, args)
  return(out)
}

#' 2D Whittaker-Henderson Smoothing
#'
#' #' Main package function to apply Whittaker-Henderson smoothing in a
#' bidimensional survival analysis framework. It takes as input a matrix of
#' observed events and a matrix of associated central exposure, both depending
#' on two covariates, and build a smooth version of the log-hazard rate.
#' Smoothing parameters may be supplied or automatically chosen according to a
#' specific criterion such as `"REML"` (the default), `"AIC"`, `"BIC"` or
#' `"GCV"`. Whittaker-Henderson may be applied in a full maximum likelihood
#' framework or an approximate gaussian framework. As Whittaker-Henderson
#' smoothing relies on full-rank smoothers, computation time and memory usage in
#' the bidimensional case may prove overwhelming and the function integrates a
#' rank-reduction procedure to avoid such issues.
#'
#' @inheritParams WH_1d
#' @param d Matrix of observed events whose rows and columns should be named.
#'   Required in case of maximum likelihood estimation
#' @param ec Matrix of central exposure. Required in case of maximum likelihood
#'   estimation. The central exposure corresponds to the sum of the exposure
#'   period over the insured population. An individual experiencing an event of
#'   interest during the year will no longer be exposed afterward and the
#'   exposure should be reduced accordingly.
#' @param lambda Smoothing parameter vector of size `2`. If missing, an
#'   optimization procedure will be used to find the optimal smoothing
#'   parameter. If provided, no optimal smoothing parameter search will take
#'   place unless the `criterion` argument is also provided, in which case
#'   `lambda` will be used as the starting parameter for the optimization
#'   procedure.
#' @param method Method to be used to find the optimal smoothing parameter.
#'   Default to `"fixed_lambda"` if `lambda` is supplied, meaning no
#'   optimization is performed. Otherwise, default to `"perf"` which means the
#'   performance iteration method is used. The alternative `"outer"` method is
#'   guaranteed to converge but slower. Both those methods rely on the `optim`
#'   function from package `stats`.
#' @param p Optional vector of size `2`. Maximum number of eigenvectors to keep
#'   on each dimension after performing the eigen decomposition of the
#'   penalization matrix. If missing, will be automatically computed so that the
#'   dimensions of (square) matrices involved in the optimization problem
#'   remains lower that the `max_dim` argument
#' @param max_dim Number of parameters to be kept in the optimization problem.
#'   Default is `200`. Values higher than `2000` may result in very high
#'   computation times and memory usage.
#' @param q Order of penalization vector of size `2`. Polynoms of degrees
#'   `(q[[1]] - 1,q[[2]] - 1)` are considered smooth and are therefore
#'   unpenalized. Should be left to the default of `c(2,2)` for most practical
#'   applications.
#' @param y Optional matrix of observations whose rows and columns should be
#'   named. Used only in the regression framework and even in this case will be
#'   automatically computed if the `d` and `ec` arguments are supplied. May be
#'   useful when using Whittaker-Henderson smoothing outside of the survival
#'   analysis framework.
#' @param wt Optional matrix of weights. As for the observation vector `y`, used
#'   only in the regression framework and even in this case will be
#'   automatically computed if the `d` argument is supplied. May be useful when
#'   using Whittaker-Henderson smoothing outside of the survival analysis
#'   framework.
#'
#' @returns An object of class `WH_2d` i.e. a list containing :
#' * `d` The inputed matrix of observed events (if supplied as input)
#' * `ec` The inputed matrix of central exposure (if supplied as input)
#' * `y` The observation matrix (either supplied or computed as y = log(d) - log(ec))
#' * `wt` The inputed matrix of weights (either supplied or computed as `d`)
#' * `y_hat` The matrix of fitted value
#' * `std_y_hat` The matrix of standard deviation associated with the fitted value
#' * `res` The matrix of model deviance residuals
#' * `edf_obs` The matrix of effective degrees of freedom associated with each observation
#' * `edf_par` The matrix of effective degrees of freedom associated with each eigenvector
#' * `diagnosis` A data.frame with one line containing the effective degrees of freedom of the model, the deviance of the fit as well as the AIC, BIC, GCV and REML criteria
#' * `Psi` The variance-covariance matrix associated with the fit, required for extrapolation.
#' * `lambda` The vector of smoothing parameters used, either supplied. or computed.
#' * `p` The number of components kept on each dimension after the rank reduction method is applied.
#' * `q` The supplied vector of orders for the penalization.
#'
#' @examples
#' keep_age <- which(rowSums(portfolio_LTC$ec) > 1e2)
#' keep_duration <- which(colSums(portfolio_LTC$ec) > 1e2)
#'
#' d  <- portfolio_LTC$d[keep_age, keep_duration]
#' ec <- portfolio_LTC$ec[keep_age, keep_duration]
#'
#' y <- log(d / ec) # observation vector
#' y[d == 0] <- - 20
#' wt <- d
#'
#' # Maximum likelihood
#' WH_2d(d, ec, lambda = c(1e2, 1e2))
#' WH_2d(d, ec) # performance iteration default method
#'
#' WH_2d(d, ec, criterion = "GCV")
#' # alternative optimization criteria for smoothing parameter selection
#'
#' # Regression
#' WH_2d(y = y, wt = wt, lambda = c(1e2, 1e2))
#' # regression framework is triggered when y is supplied
#' WH_2d(d, y = y, lambda = c(1e2, 1e2)) # d is used as wt
#' WH_2d(d, ec, framework = "reg", lambda = c(1e2, 1e2))
#' # setting framework = "reg" forces computation of y from d and ec
#'
#' keep_age2 <- which(rowSums(portfolio_LTC$ec) > 0)
#' keep_duration2 <- which(colSums(portfolio_LTC$ec) > 0)
#'
#' # Rank reduction
#' d  <- portfolio_LTC$d[keep_age2, keep_duration2]
#' ec <- portfolio_LTC$ec[keep_age2, keep_duration2]
#'
#' prod(dim(d)) # problem dimension is 1,232 !
#' WH_2d(d, ec)
#' # rank-reduction is used to find an approximate solution with 200 parameters
#'
#' @export
WH_2d <- function(d, ec, lambda, criterion, method, max_dim = 200, p,
                  q = c(2, 2), framework, y, wt, quiet = FALSE, ...) {

  if (missing(framework)) framework <- if (!missing(y)) "reg" else "ml"
  if (length(framework) != 1 || !(framework %in% c("reg", "ml"))) {
    stop("framework should be one of ml (the default) or reg")
  }
  if (framework == "reg") {
    if (missing(y)) {
      if (!missing(d) && !missing(ec)) {
        if (is.null(dim(d)) || length(dim(d)) != 2) stop("d should be a numeric matrix")
        if (is.null(dim(ec)) || length(dim(ec)) != 2) stop("ec should be a numeric matrix")
        if (any(dim(d) != dim(ec))) stop("Dimensions of d and ec must match")
        if (!quiet) message("Computing y as y = log(d / ec)")
        y <- log(d / ec)
        y[d == 0] <- - 20
      } else {
        stop("Either y or both d and ec required in the regression framework")
      }
    } else {
      if (is.null(dim(y)) || length(dim(y)) != 2) stop("y should be a numeric matrix")
    }
    if (missing(wt)) {
      if (!missing(d)) {
        if (!quiet) message("Using d as weights")
        wt <- d
      } else {
        stop("Either wt or d needs to be supplied in the regression framework")
      }
    } else {
      if (is.null(dim(wt)) || length(dim(wt)) != 2) stop("wt should be a numeric matrix")
    }
    if (any(dim(y) != dim(wt))) stop("Dimensions of y and wt must match")
    if (is.null(rownames(y))) {
      y_rownames <- seq_along(dim(y)[[1]])
      if (!quiet) warning("No names found for y rows, setting names to: ", paste(range(y_rownames), collapse = " - "))
    }
    if (is.null(colnames(y))) {
      y_colnames <- seq_along(dim(y)[[2]])
      if (!quiet) warning("No names found for y columns, setting names to: ", paste(range(y_colnames), collapse = " - "))
    }
  } else {
    if (missing(d)) stop("d argument required in the maximum likelihood framework")
    if (missing(ec)) stop("ec argument required in the maximum likelihood framework")
    if (is.null(dim(d)) || length(dim(d)) != 2) stop("d should be a numeric matrix")
    if (is.null(dim(ec)) || length(dim(ec)) != 2) stop("ec should be a numeric matrix")
    if (any(dim(d) != dim(ec))) stop("Dimensions of d and ec must match")
    if (is.null(rownames(d))) {
      d_rownames <- seq_along(dim(d)[[1]])
      if (!quiet) warning("No names found for d rows, setting names to: ", paste(range(d_rownames), collapse = " - "))
    }
    if (is.null(colnames(d))) {
      d_colnames <- seq_along(dim(d)[[2]])
      if (!quiet) warning("No names found for d columns, setting names to: ", paste(range(d_colnames), collapse = " - "))
    }
  }

  if (missing(criterion)) criterion <- "REML"
  criteria <- c("REML", "AIC", "BIC", "GCV")
  if (length(criterion) != 1 || !(criterion %in% criteria)) stop(
    "criterion should be one of ", paste(criteria, collapse = ", "))

  if (missing(method)) {
    if (!missing(lambda)) {
      if (length(lambda) != 2) stop("smoothing parameter vector lambda should be of length 2")
      if (!is.numeric(lambda) || any(lambda <= 0)) stop(
        "smoothing parameters should be strictly positive numbers")
      method <- "fixed_lambda"
      max_dim <- Inf
    } else {
      if (!quiet) message("Using performance iteration / Nelder-Mead method")
      method <- "perf"
    }
  }
  methods <- c("fixed_lambda", "perf", "outer")
  if (length("method") != 1 || !(method %in% methods)) stop(
    "method should be one of ", paste(methods, collapse = ", "))

  if (missing(lambda)) {
    lambda <- c(1e3, 1e3)
  } else {
    if (method != "fixed_lambda" && !quiet) message(
      "Both method and lambda specified, lambda used as starting parameter")
  }

  n <- if (framework == "reg") dim(y) else dim(d)
  if (missing(p)) {

    max_ratio <- sqrt(max_dim / (n[[2]] * n[[1]]))
    p <- floor(pmin(max_ratio, 1) * n)
  } else {
    if (!is.numeric(q) || length(q) != 2 || any(q <= 0) ||
        (max(abs(q - round(q))) > .Machine$double.eps ^ 0.5)) stop(
          "p should be an integer vector of length 2")
    if (any(p <= q) || any(p > n)) stop(
      "p should range between ", paste(q, collapse = " / "), " and ", paste(n, collapse = " / "))
  }

  if (!is.numeric(q) || length(q) != 2 || any(q <= 0) ||
      (max(abs(q - round(q))) > .Machine$double.eps ^ 0.5)) stop(
        "q should be a couple of strictly positive integers"
      )
  if (any(q > 3) && !quiet) warning("Differences of order q > 3 have no meaningful interpretation.\n",
                                    "Consider using a lower value for q")

  reg <- (framework == "reg")
  what <- paste("WH_2d", method, sep = "_")
  args <- (if (reg) list(y = y, wt = wt) else list(d = d, ec = ec)) |>
    c(list(p = p, q = q, lambda = lambda, reg = reg),
      if (method %in% c("perf", "outer")) list(criterion = criterion),
      list(...))

  out <- do.call(what, args)
  return(out)
}

# Extrapolation----

#' Predict new values using a fitted 1D WH model
#'
#' Extrapolate the model for new observations.
#'
#' @param object An object of class `"WH_1d"` returned by the [WH_1d()] function
#' @param newdata A vector containing the position of new observations.
#'   Observations from the fit will automatically be added to this, in the
#'   adequate order
#' @param ... Not used
#'
#' @returns An object of class `"WH_1d"` with additional components for model
#'   prediction.
#'
#' @examples
#' keep <- which(portfolio_mort$ec > 0)
#' d <- portfolio_mort$d[keep]
#' ec <- portfolio_mort$ec[keep]
#'
#' WH_1d(d, ec) |> predict(newdata = 18:99) |> plot()
#'
#' @export
predict.WH_1d <- function(object, newdata = NULL, ...) {

  if (!inherits(object, "WH_1d")) stop("object must be of class WH_1d")
  if (!is.numeric(newdata)) stop("newdata should be a vector containing the names of predicted values")

  data <- as.numeric(names(object$y))
  full_data <- sort(union(data, newdata))

  n <- length(data)
  n_pred <- length(full_data)
  n_new <- n_pred - n

  ind_fit <- which(full_data %in% data)
  ind_inf <- which(full_data < min(data))
  ind_sup <- which(full_data > max(data))
  ind_new <- c(ind_inf, ind_sup)

  D_mat_pred <- build_D_mat(n_pred, object$q)

  D1_inf <- D_mat_pred[ind_inf, ind_fit, drop = FALSE]
  D1_sup <- D_mat_pred[ind_sup - object$q, ind_fit, drop = FALSE]
  D1 <- rbind(D1_inf, D1_sup)

  D2_inf <- D_mat_pred[ind_inf, ind_inf, drop = FALSE]
  D2_sup <- D_mat_pred[ind_sup - object$q, ind_sup, drop = FALSE]
  D2 <- blockdiag(D2_inf, D2_sup)

  D2_inv <- if (n_new == 0) matrix(0, 0, 0) else solve(D2)

  A_pred <- matrix(0, n_pred, n)
  A_pred[ind_fit,] <- object$U
  A_pred[ind_new,] <- - D2_inv %*% D1 %*% object$U

  y_pred <- c(A_pred %*% object$beta_hat)
  std_y_pred <- sqrt(rowSums(A_pred * (A_pred %*% object$Psi)))

  names(y_pred) <- names(std_y_pred) <- full_data

  object$y_pred <- y_pred
  object$std_y_pred <- std_y_pred

  return(object)
}

#' Predict new values using a fitted 2D WH model
#'
#' Extrapolate the model for new observations in a way that is consistent with
#' the fitted values
#'
#' @param object An object of class `"WH_2d"` returned by the [WH_2d()] function
#' @param newdata A list containing two vectors indicating the new observation
#'   positions
#' @param ... Not used
#'
#' @returns An object of class `"WH_2d"` with additional components for model
#'   prediction.
#'
#' @examples
#' keep_age <- which(rowSums(portfolio_LTC$ec) > 1e2)
#' keep_duration <- which(colSums(portfolio_LTC$ec) > 1e2)
#'
#' d  <- portfolio_LTC$d[keep_age, keep_duration]
#' ec <- portfolio_LTC$ec[keep_age, keep_duration]
#'
#' WH_2d(d, ec) |> predict(newdata = list(age = 50:99, duration = 0:19)) |> plot()
#'
#' @export
predict.WH_2d <- function(object, newdata = NULL, ...) {

  if (!inherits(object, "WH_2d")) stop("object must be of class WH_2d")
  if (length(newdata) != 2 || !is.numeric(newdata[[1]]) || !is.numeric(newdata[[2]])) stop(
    "newdata should be a list with two elements containing the row names and column names for predicted values")

  data <- dimnames(object$y) |> purrr::map(as.numeric)
  full_data <- purrr::map2(data, newdata, \(x,y) sort(union(x, y)))
  ind_fit <- purrr::map2(data, full_data, \(x,y) which(y %in% x))

  n <- purrr::map_int(data, length)
  n_inf <- purrr::map2_int(data, full_data, \(x,y) sum(y < min(x)))
  n_sup <- purrr::map2_int(data, full_data, \(x,y) sum(y > max(x)))
  n_pred <- n + n_inf + n_sup

  wt_pred <- matrix(0, n_pred[[1]], n_pred[[2]])
  wt_pred[ind_fit[[1]], ind_fit[[2]]] <- object$wt

  # W_pred <- diag(wt_pred) # extended weight matrix
  D_mat_pred <- purrr::map2(n_pred, object$q, build_D_mat) # extended difference matrices
  P_pred <- object$lambda[[1]] * diag(n_pred[[2]]) %x% crossprod(D_mat_pred[[1]]) +
    object$lambda[[2]] * crossprod(D_mat_pred[[2]]) %x% diag(n_pred[[1]]) # extended penalization matrix
  diag(P_pred) <- diag(P_pred) + c(wt_pred)
  Psi_pred <- P_pred |> chol() |> chol2inv() # unconstrained variance / covariance matrix

  ind_rows <- c(rep(FALSE, n_inf[[1]]), rep(TRUE, n[[1]]), rep(FALSE, n_sup[[1]]))
  ind_coef_2d <- c(rep(FALSE, n_pred[[1]] * n_inf[[2]]), rep(ind_rows, n[[2]])) |> which()

  Psi_inv <- Psi_pred[ind_coef_2d, ind_coef_2d] |> chol() |> chol2inv()
  A_pred <- Psi_pred[,ind_coef_2d] %*% Psi_inv

  y_pred <- c(A_pred %*% c(object$y_hat))
  std_y_pred <- sqrt(rowSums(A_pred * (A_pred %*% (object$U %*% object$Psi %*% t(object$U)))))

  dim(y_pred) <- dim(std_y_pred) <- n_pred # set dimension for output matrices
  dimnames(y_pred) <- dimnames(std_y_pred) <- full_data # set names for output matrices

  object$y_pred <- y_pred
  object$std_y_pred <- std_y_pred

  return(object)
}

# Formatting----

#' Store WH model fit results in a data.frame
#'
#' @param object An object of class  `"WH_1d"` or `"WH_2d"` returned
#'   by one of the eponymous functions [WH_1d()] or [WH_2d()]
#' @param dim1 The (optional) name to be given to the first dimension
#' @param dim2 The (optional) name to be given to the second dimension
#'
#' @returns A data.frame gathering information about the fitted and predicted
#'   values, the model variance, residuals and effective degrees of freedom...
#'
#' @export
output_to_df <- function(object, dim1 = "x", dim2 = "t") {

  if (!inherits(object, c("WH_1d", "WH_2d"))) stop("object must be of class WH_1d or WH_2d")
  if (length(dim1) != 1 || length(dim2) != 1) stop("The optional dim1 and dim2 arguments should be of length 1")

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

    df$y <- df$wt <- df$res <- df$edf_obs <- NA_real_
    df$y[x_pred %in% x] <- c(object$y)
    df$wt[x_pred %in% x] <- c(object$wt)
    df$res[x_pred %in% x] <- c(object$res)
    df$edf_obs[x_pred %in% x] <- c(object$edf_obs)

    if("ec" %in% names(object) && "d" %in% names(object)) {

      df$ec <- df$d <- NA_real_

      df$ec[x_pred %in% x] <- c(object$ec)
      df$d[x_pred %in% x] <- c(object$d)
    }

    df[,c(dim1, intersect(
      c("ec", "d", "y", "y_hat", "std_y_hat", "wt", "res", "edf_obs"),
      names(object)))]

  } else {

    x <- as.numeric(rownames(object$y))
    x_pred <- as.numeric(rownames(object$y_hat))

    t <- as.numeric(colnames(object$y))
    t_pred <- as.numeric(colnames(object$y_hat))

    df[[dim1]] <- rep(x_pred, times = length(t_pred))
    df[[dim2]] <- rep(t_pred, each = length(x_pred))

    df$y <- df$wt <- df$res <- df$edf_obs <- NA_real_
    df$y[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$y)
    df$wt[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$wt)
    df$res[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$res)
    df$edf_obs[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$edf_obs)

    if("ec" %in% names(object) && "d" %in% names(object)) {

      df$ec <- df$d <- NA_real_

      df$ec[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$ec)
      df$d[df[[dim1]] %in% x & df[[dim2]] %in% t] <- c(object$d)
    }

    df[,c(dim1, dim2, intersect(
      c("ec", "d", "y", "y_hat", "std_y_hat", "wt", "res", "edf_obs"),
      names(object)))]
  }
}

# Display----

#' Display of 1D WH object
#'
#' @param x An object of class `"WH_1d"` returned by the [WH_1d()] function
#' @param ... Not used
#'
#' @returns Invisibly returns `x`.
#'
#' @examples
#' keep <- which(portfolio_mort$ec > 0)
#' d <- portfolio_mort$d[keep]
#' ec <- portfolio_mort$ec[keep]
#'
#' y <- log(d / ec)
#' y[d == 0] <- - 20
#' wt <- d
#'
#' WH_1d(d, ec) |> summary()
#' WH_1d(d, ec, method = "perf") |> summary()
#' WH_1d(d, ec, criterion = "GCV") |> summary()
#' WH_1d(y = y, wt = wt) |> summary()
#'
#' @export
print.WH_1d <- function(x, ...) {

  if (!inherits(x, "WH_1d")) stop("x must be of class WH_2d")

  cat("An object fitted using the WH_1D function\n")
  cat("Initial data contains", length(x$y), "data points:\n")
  cat("  Observation positions: ", as.numeric(names(x$y)[[1]]), " to ",as.numeric(names(x$y)[[length(x$y)]]), "\n")
  cat("Smoothing parameter selected:", format(x$lambda, digits = 2), "\n")
  cat("Associated degrees of freedom:", format(sum(x$edf_obs), digits = 2), "\n\n")
  invisible(x)
}

#' Display of 2D WH object
#'
#' @param x An object of class `"WH_2d"` returned by the [WH_2d()] function
#' @param ... Not used
#'
#' @returns Invisibly returns `x`.
#'
#' @examples
#' keep_age <- which(rowSums(portfolio_LTC$ec) > 1e2)
#' keep_duration <- which(colSums(portfolio_LTC$ec) > 1e2)
#'
#' d  <- portfolio_LTC$d[keep_age, keep_duration]
#' ec <- portfolio_LTC$ec[keep_age, keep_duration]
#'
#' WH_2d(d, ec) |> summary()
#'
#' @export
print.WH_2d <- function(x, ...) {

  if (!inherits(x, "WH_2d")) stop("x must be of class WH_2d")

  cat("An object fitted using the WH_2D function\n")
  cat("Initial data contains", prod(dim(x$y)), "data points:\n")
  cat("  First dimension: ", as.numeric(rownames(x$y)[[1]]), " to ",as.numeric(rownames(x$y)[[nrow(x$y)]]), "\n")
  cat("  Second dimension: ", as.numeric(colnames(x$y)[[1]]), " to ",as.numeric(colnames(x$y)[[ncol(x$y)]]), "\n")
  cat("Smoothing parameters selected:", format(x$lambda, digits = 2), "\n")
  cat("Associated degrees of freedom:", format(sum(x$edf_obs), digits = 2), "\n\n")
  invisible(x)
}

# Plots----

#' Plot 1D WH fit
#'
#' @param x An object of class `"WH_1d"` returned by the [WH_1d()] function
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
#' keep <- which(portfolio_mort$ec > 0)
#' d <- portfolio_mort$d[keep]
#' ec <- portfolio_mort$ec[keep]
#'
#' WH_1d(d, ec) |> plot()
#' WH_1d(d, ec) |> plot("res")
#' WH_1d(d, ec) |> plot("edf")
#'
#' @export
plot.WH_1d <- function(x, what = "fit", trans, ...) {

  if (!inherits(x, "WH_1d")) stop("x must be of class WH_2d")
  whats <- c("fit", "res", "edf")
  if (length(what) != 1 || !(what %in% whats)) stop(
    "what should be one of ", paste(whats, collapse = ", "))

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
           plot(df$x, df$edf_obs,
                xlab = "x", ylab = "degrees of freedom", type = "b")
           graphics::abline(a = 0, b = 0, lty = 2, col = "blue")
           graphics::abline(a = 1, b = 0, lty = 2, col = "blue")
         })
}

#' Plot 2D WH fit
#'
#' @param x An object of class `"WH_2d"` returned by the [WH_2d()] function
#' @param what What should be plotted (y_hat, std_y_hat, res, edf)
#' @param trans An (optional) transformation to be applied to the data
#' @param ... Not used
#'
#' @returns A plot representing the given element of the fit...
#'
#' @examples
#'
#' keep_age <- which(rowSums(portfolio_LTC$ec) > 1e2)
#' keep_duration <- which(colSums(portfolio_LTC$ec) > 1e2)
#'
#' d  <- portfolio_LTC$d[keep_age, keep_duration]
#' ec <- portfolio_LTC$ec[keep_age, keep_duration]
#'
#' WH_2d(d, ec) |> plot()
#' WH_2d(d, ec) |> plot("std_y_hat")
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

# 1D Fit----

#' Whittaker-Henderson Smoothing (Maximum Likelihood, fixed lambda)
#'
#' @param d Vector of observed events
#' @param ec Vector of central exposure
#' @param y Vector of observations
#' @param wt Optional vector of weights
#' @param lambda Smoothing parameter
#' @param p The number of eigenvectors to keep
#' @param q Order of penalization. Polynoms of degrees q - 1 are considered
#'   smooth and are therefore unpenalized
#' @param reg Should the regression framework be used ? Boolean. If `TRUE`, will
#'   stop after the first iteration.
#' @param accu_dev Tolerance for the convergence of the optimization procedure
#'
#' @returns An object of class `"WH_1d"` i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#' @keywords internal
WH_1d_fixed_lambda <- function(d, ec, y, wt, lambda = 1e3, q = 2, p,
                               reg = FALSE, verbose = FALSE, accu_dev = 1e-12) {

  # Initialization
  n <- if (reg) length(y) else length(d)
  if (missing(p)) p <- n
  eig <- eigen_dec(n, q, p)
  U <- eig$U
  s <- eig$s

  sum_wt <- if (reg) sum(wt) else sum(d)
  which_pos <- if (reg) which(wt != 0) else which(ec != 0)
  n_pos <- length(which_pos)
  U_pos <- U[which_pos,]

  if (reg) {

    z <- y
    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]
    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

  } else {

    off <- log(pmax(ec, 1e-4))
    y <- ifelse(d == 0, NA, log(d)) - off
    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
  }

  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    old_dev_pen <- dev_pen

    # update of parameters, working vector and weight matrix
    if (!reg) {

      wt <- new_wt
      z <- y_hat + d / wt - 1
      wt_pos <- c(wt)[which_pos]
      z_pos <- c(z)[which_pos]
      tUWU <- t(U_pos) %*% (wt_pos * U_pos)
      tUWz <- t(U_pos) %*% (wt_pos * z_pos)
    }

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + lambda * s
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    beta_hat <- c(Psi %*% tUWz) # fitted parameter

    y_hat <- c(U %*% beta_hat)
    if (!reg) new_wt <- exp(y_hat + off)
    res <- if (reg) (sqrt(wt) * (y - y_hat)) else compute_res_deviance(d, new_wt) # (weighted) residuals
    dev <- sum(res * res)

    RESS <- sum(beta_hat * s * beta_hat)
    pen <- lambda * RESS

    dev_pen <- dev + pen

    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- if (reg) FALSE else ((old_dev_pen - dev_pen) > accu_dev * sum_wt)
  }

  edf_par <- colSums(Psi * tUWU) # effective degrees of freedom by parameter

  aux <- rowSums(U * (U %*% Psi))
  edf_obs <- wt * aux # effective degrees of freedom by observation / parameter
  std_y_hat <- sqrt(aux) # standard deviation of fit

  sum_edf <- sum(edf_par) # effective degrees of freedom
  tr_log_P <- (p - q) * log(lambda) + sum(log(s[- seq_len(q)]))
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  diagnosis <- get_diagnosis(dev, pen, sum_edf, n_pos, tr_log_P, tr_log_Psi)

  names(y_hat) <- names(std_y_hat) <- names(res) <- names(edf_obs) <-
    names(wt) <- names(z) <- names(y) # set names for output vectors

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf_obs = edf_obs,
              beta_hat = beta_hat, edf_par = edf_par, diagnosis = diagnosis,
              U = U, Psi = Psi, lambda = lambda, p = p, q = q)
  class(out) <- "WH_1d"

  return(out)
}

#' Whittaker-Henderson Smoothing (Maximum Likelihood, optimize function)
#'
#' @inheritParams WH_1d_fixed_lambda
#' @param criterion Criterion used to choose the smoothing parameter. One of
#'   "GCV" (default), "AIC" or "BIC".
#' @param lambda Initial smoothing parameter
#' @param verbose Should information about the optimization progress be
#'   displayed
#' @param accu_crit Tolerance for the convergence of the outer optimization
#'   procedure
#'
#' @returns An object of class `"WH_1d"` i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#' @keywords internal
WH_1d_outer <- function(d, ec, y, wt, q = 2, p, criterion = "REML", lambda = 1e3,
                        reg = FALSE, verbose = FALSE, accu_crit = 1e-12, accu_dev = 1e-12) {

  # Initialization
  n <- if (reg) length(y) else length(d)
  if (missing(p)) p <- n
  eig <- eigen_dec(n, q, p)
  U <- eig$U
  s <- eig$s

  sum_wt <- if (reg) sum(wt) else sum(d)
  which_pos <- if (reg) which(wt != 0) else which(ec != 0)
  n_pos <- length(which_pos)
  U_pos <- U[which_pos,]

  if (reg) {

    z <- y
    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]
    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

  } else {

    off <- log(pmax(ec, 1e-4))
    y <- ifelse(d == 0, NA, log(d)) - off
    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
  }

  WH_1d_aux <- function(log_lambda) {

    lambda <- exp(log_lambda)

    if (!reg) {

      y_hat <- log(pmax(d, 1e-8)) - off
      new_wt <- exp(y_hat + off)
    }

    dev_pen <- Inf
    cond_dev_pen <- TRUE

    # Loop
    while (cond_dev_pen) {

      old_dev_pen <- dev_pen

      # update of working vector and weight matrix
      if (!reg) {

        wt <- new_wt
        z <- y_hat + d / wt - 1
        wt_pos <- c(wt)[which_pos]
        z_pos <- c(z)[which_pos]
        tUWU <- t(U_pos) %*% (wt_pos * U_pos)
        tUWz <- t(U_pos) %*% (wt_pos * z_pos)
      }

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + lambda * s
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      beta_hat <- c(Psi %*% tUWz) # fitted value

      y_hat <- c(U %*% beta_hat)
      if (!reg) new_wt <- exp(y_hat + off)
      res <- if (reg) (sqrt(wt) * (y - y_hat)) else compute_res_deviance(d, new_wt) # (weighted) residuals
      dev <- sum(res * res)

      RESS <- sum(beta_hat * s * beta_hat)
      pen <- lambda * RESS

      dev_pen <- dev + pen

      if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                       "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
      cond_dev_pen <-  if (reg) FALSE else (old_dev_pen - dev_pen) > accu_dev * sum_wt
    }

    if (criterion != "REML") sum_edf <- sum(Psi * tUWU) # effective degrees of freedom

    score <- switch(criterion,
                    AIC = dev + 2 * sum_edf,
                    BIC = dev + log(n_pos) * sum_edf,
                    GCV = n_pos * dev / (n_pos - sum_edf) ^ 2,
                    REML = {
                      tr_log_P <- (p - q) * log(lambda) + sum(log(s[- seq_len(q)]))
                      tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())
                      REML <- dev + pen - tr_log_P + tr_log_Psi
                    })

    if (verbose) cat(criterion, " :", format(
      if (criterion == "REML") - score / 2 else score, digits = 3, decimal.mark = ","), "\n")
    return(score)
  }

  lambda <- exp(stats::optimize(f = WH_1d_aux, interval = 25 * c(- 1, 1), tol = accu_crit * sum_wt)$minimum)
  out <- WH_1d_fixed_lambda(d = d, ec = ec, y = y, wt = wt, lambda = lambda, p = n, q = q, reg = reg,
                            verbose = verbose, accu_dev = accu_dev)

  return(out)
}

#' Whittaker-Henderson Smoothing (Maximum Likelihood, Generalized Fellner-Schall update)
#'
#' @inheritParams WH_1d_outer
#'
#' @returns An object of class `"WH_1d"` i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#' @keywords internal
WH_1d_perf <- function(d, ec, y, wt, q = 2, p, criterion = "REML", lambda = 1e3,
                       reg = FALSE, verbose = FALSE, accu_crit = 1e-12, accu_dev = 1e-12) {

  # Initialization
  n <- if (reg) length(y) else length(d)
  if (missing(p)) p <- n
  eig <- eigen_dec(n, q, p)
  U <- eig$U
  s <- eig$s

  sum_wt <- if (reg) sum(wt) else sum(d)
  which_pos <- if (reg) which(wt != 0) else which(ec != 0)
  n_pos <- length(which_pos)
  U_pos <- U[which_pos,]

  if (reg) {

    z <- y
    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]
    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

  } else {

    off <- log(pmax(ec, 1e-4))
    y <- ifelse(d == 0, NA, log(d)) - off
    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
  }

  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    old_dev_pen <- dev_pen

    # update of working vector and weight matrix
    if (!reg) {

      wt <- new_wt
      z <- y_hat + d / wt - 1
      wt_pos <- c(wt)[which_pos]
      z_pos <- c(z)[which_pos]
      tUWU <- t(U_pos) %*% (wt_pos * U_pos)
      tUWz <- t(U_pos) %*% (wt_pos * z_pos)
    }

    WH_1d_aux <- function(log_lambda) {

      lambda <- exp(log_lambda)

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + lambda * s
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      beta_hat <- c(Psi %*% tUWz) # fitted value

      y_hat <- c(U %*% beta_hat)
      res <- sqrt(wt) * ((if (reg) y else z) - y_hat) # (weighted) residuals
      dev <- sum(res * res)

      if (criterion != "REML") sum_edf <- sum(Psi * tUWU) # effective degrees of freedom

      score <- switch(criterion,
                      AIC = dev + 2 * sum_edf,
                      BIC = dev + log(n_pos) * sum_edf,
                      GCV = n_pos * dev / (n_pos - sum_edf) ^ 2,
                      REML = {
                        RESS <- sum(beta_hat * s * beta_hat)
                        pen <- lambda * RESS
                        tr_log_P <- (p - q) * log(lambda) + sum(log(s[- seq_len(q)]))
                        tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())
                        REML <- dev + pen - tr_log_P + tr_log_Psi
                      })

      if (verbose) cat(criterion, " :", format(
        if (criterion == "REML") - score / 2 else score, digits = 3, decimal.mark = ","), "\n")

      return(score)
    }

    lambda <- exp(stats::optimize(f = WH_1d_aux, interval = 25 * c(- 1, 1), tol = accu_crit * sum_wt)$minimum)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + lambda * s
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    beta_hat <- c(Psi %*% tUWz) # fitted value

    y_hat <- c(U %*% beta_hat)
    if (!reg) new_wt <- exp(y_hat + off)
    res <- if (reg) (sqrt(wt) * (y - y_hat)) else compute_res_deviance(d, new_wt) # (weighted) residuals
    dev <- sum(res * res)

    RESS <- sum(beta_hat * s * beta_hat)
    pen <- lambda * RESS

    dev_pen <- dev + pen

    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <-  if (reg) FALSE else (old_dev_pen - dev_pen) > accu_dev * sum_wt
  }

  out <- WH_1d_fixed_lambda(d = d, ec = ec, y = y, wt = wt, lambda = lambda, p = n, q = q, reg = reg,
                            verbose = verbose, accu_dev = accu_dev)

  return(out)
}

# 2D Fit----

#' 2D Whittaker-Henderson Smoothing (Maximum Likelihood, fixed lambda)
#'
#' @param y Matrix of observations
#' @param wt Optional matrix of weights
#' @param lambda Vector of smoothing parameter
#' @param p The number of eigenvectors to keep on each dimension
#' @param q Matrix of orders of penalization. Polynoms of degrees q - 1 are considered
#'   smooth and are therefore unpenalized
#'
#' @returns An object of class `"WH_2d"` i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#' @keywords internal
WH_2d_fixed_lambda <- function(d, ec, y, wt, lambda = c(1e3, 1e3), q = c(2, 2), p,
                               reg = FALSE, verbose = FALSE, accu_dev = 1e-12) {

  # Initialization
  n <- if (reg) dim(y) else dim(d)
  if (missing(p)) p <- n
  eig <- purrr::pmap(list(n = n, q = q, p = p), eigen_dec)
  U_eig <- purrr::map(eig, "U")
  U <- U_eig |> rev() |> purrr::reduce(kronecker)
  s_eig <- purrr::map(eig, "s")
  s <- list(rep(s_eig[[1]], p[[2]]), rep(s_eig[[2]], each = p[[1]]))

  sum_wt <- if (reg) sum(wt) else sum(d)
  which_pos <- if (reg) which(wt != 0) else which(ec != 0)
  n_pos <- length(which_pos)
  U_pos <- U[which_pos,]

  if (reg) {

    z <- y
    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]
    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

  } else {

    off <- log(pmax(ec, 1e-4))
    y <- ifelse(d == 0, NA, log(d)) - off
    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
  }

  s_lambda <- purrr::map2(lambda, s, `*`)
  sum_s_lambda <- s_lambda |> do.call(what = `+`)

  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    old_dev_pen <- dev_pen

    # update of parameters, working vector and weight matrix
    if (!reg) {

      wt <- new_wt
      z <- y_hat + d / wt - 1
      wt_pos <- c(wt)[which_pos]
      z_pos <- c(z)[which_pos]
      tUWU <- t(U_pos) %*% (wt_pos * U_pos)
      tUWz <- t(U_pos) %*% (wt_pos * z_pos)
    }

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    beta_hat <- c(Psi %*% tUWz) # fitted value
    y_hat <- c(U %*% beta_hat)
    if (!reg) new_wt <- exp(y_hat + off)
    res <- if (reg) sqrt(wt) * (y - y_hat) else compute_res_deviance(d, new_wt) # (weighted) residuals
    dev <- sum(res * res)

    RESS <- purrr::map_dbl(s, \(x) sum(beta_hat * x * beta_hat))
    pen <- purrr::map2(lambda, RESS, `*`) |> do.call(what = `+`)

    dev_pen <- dev + pen

    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- if (reg) FALSE else (old_dev_pen - dev_pen) > accu_dev * sum_wt
  }

  edf_par <- colSums(Psi * tUWU) |> matrix(p[[1]], p[[2]]) # effective degrees of freedom by parameter
  omega_j <- purrr::map(s_lambda, \(x) ifelse(x == 0, 0, x / sum_s_lambda))
  sum_edf_random <- purrr::map_dbl(omega_j, \(x) sum(x * edf_par))

  aux <- rowSums(U * (U %*% Psi))
  edf_obs <-c(wt) * aux # effective degrees of freedom by observation / parameter
  std_y_hat <- sqrt(aux) # standard deviation of fit

  sum_edf <- sum(edf_par) # effective degrees of freedom
  tr_log_P <- purrr::map2(lambda, s, `*`) |> do.call(what = `+`) |> Filter(f = \(x) x > 0) |> log() |> sum()
  tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())

  diagnosis <- get_diagnosis(dev, pen, sum_edf, n_pos, tr_log_P, tr_log_Psi)

  dim(y_hat) <- dim(std_y_hat) <- dim(res) <- dim(edf_obs) <-
    dim(wt) <- dim(y) # set dimensions for output matrices
  dimnames(y_hat) <- dimnames(std_y_hat) <- dimnames(res) <- dimnames(edf_obs) <-
    dimnames(wt) <- dimnames(y) # set names for output matrices

  out <- list(y = y, wt = wt, y_hat = y_hat, std_y_hat = std_y_hat, res = res, edf_obs = edf_obs,
              beta_hat = beta_hat, edf_par = edf_par, diagnosis = diagnosis,
              U = U, Psi = Psi, lambda = lambda, p = p, q = q)
  class(out) <- "WH_2d"

  return(out)
}

#' 2D Whittaker-Henderson Smoothing (Maximum Likelihood, optim function)
#'
#' @inheritParams WH_1d_outer
#' @inheritParams WH_2d_fixed_lambda
#' @param lambda Initial smoothing parameters
#'
#' @returns An object of class `"WH_2d"` i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#' @keywords internal
WH_2d_outer <- function(d, ec, y, wt, q = c(2, 2), p, criterion = "REML", lambda = c(1e3, 1e3),
                        reg = FALSE, verbose = FALSE, accu_crit = 1e-12, accu_dev = 1e-12) {

  # Initialization
  n <- if (reg) dim(y) else dim(d)
  if (missing(p)) p <- n
  eig <- purrr::pmap(list(n = n, q = q, p = p), eigen_dec)
  U_eig <- purrr::map(eig, "U")
  U <- U_eig |> rev() |> purrr::reduce(kronecker)
  s_eig <- purrr::map(eig, "s")
  s <- list(rep(s_eig[[1]], p[[2]]), rep(s_eig[[2]], each = p[[1]]))

  sum_wt <- if (reg) sum(wt) else sum(d)
  which_pos <- if (reg) which(wt != 0) else which(ec != 0)
  n_pos <- length(which_pos)
  U_pos <- U[which_pos,]

  if (reg) {

    z <- y
    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]
    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

  } else {

    off <- log(pmax(ec, 1e-4))
    y <- ifelse(d == 0, NA, log(d)) - off
    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
  }

  WH_2d_aux <- function(log_lambda) {

    lambda <- exp(log_lambda)
    s_lambda <- purrr::map2(lambda, s, `*`)
    sum_s_lambda <- s_lambda |> do.call(what = `+`)

    if (!reg) {

      y_hat <- log(pmax(d, 1e-8)) - off
      new_wt <- exp(y_hat + off)
    }

    dev_pen <- Inf
    cond_dev_pen <- TRUE

    # Loop
    while (cond_dev_pen) {

      old_dev_pen <- dev_pen

      # update of parameters, working vector and weight matrix
      if (!reg) {

        wt <- new_wt
        z <- y_hat + d / wt - 1
        wt_pos <- c(wt)[which_pos]
        z_pos <- c(z)[which_pos]
        tUWU <- t(U_pos) %*% (wt_pos * U_pos)
        tUWz <- t(U_pos) %*% (wt_pos * z_pos)
      }

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      beta_hat <- c(Psi %*% tUWz) # fitted value
      y_hat <- c(U %*% beta_hat)
      if (!reg) new_wt <- exp(y_hat + off)

      res <- if (reg) sqrt(wt) * (y - y_hat) else compute_res_deviance(d, new_wt) # (weighted) residuals
      dev <- sum(res * res)

      RESS <- purrr::map_dbl(s, \(x) sum(beta_hat * x * beta_hat))
      pen <- purrr::map2(lambda, RESS, `*`) |> do.call(what = `+`)

      dev_pen <- dev + pen

      if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                       "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
      cond_dev_pen <- if (reg) FALSE else (old_dev_pen - dev_pen) > accu_dev * sum_wt
    }

    sum_edf <- if (criterion != "REML") sum(Psi * tUWU) # effective degrees of freedom

    score <- switch(criterion,
                    AIC = dev + 2 * sum_edf,
                    BIC = dev + log(prod(n_pos)) * sum_edf,
                    GCV = prod(n_pos) * dev / (prod(n_pos) - sum_edf) ^ 2,
                    REML = {
                      tr_log_P <- purrr::map2(lambda, s, `*`) |> do.call(what = `+`) |> Filter(f = \(x) x > 0) |> log() |> sum()
                      tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())
                      REML <- dev + pen - tr_log_P + tr_log_Psi
                    })

    if (verbose) cat(criterion, " :", format(
      if (criterion == "REML") - score / 2 else score, digits = 3, decimal.mark = ","), "\n")

    return(score)
  }

  lambda <- exp(stats::optim(par = log(lambda),
                             fn = WH_2d_aux,
                             control = list(reltol = accu_crit * sum_wt))$par)
  out <- WH_2d_fixed_lambda(d = d, ec = ec, y = y, wt = wt, lambda = lambda, p = n, q = q, reg = reg,
                            verbose = verbose, accu_dev = accu_dev)

  return(out)
}

#' 2D Whittaker-Henderson Smoothing (Maximum Likelihood, Generalized Fellner-Schall update)
#'
#' @inheritParams WH_2d_outer
#'
#' @returns An object of class `"WH_2d"` i.e. a list containing model fit,
#'   variance, residuals and degrees of freedom as well as diagnosis to asses
#'   the quality of the fit.
#' @keywords internal
WH_2d_perf <- function(d, ec, y, wt, q = c(2, 2), p, criterion = "REML", lambda = c(1e3, 1e3),
                       reg = FALSE, verbose = FALSE, accu_crit = 1e-12, accu_dev = 1e-12) {

  # Initialization
  n <- if (reg) dim(y) else dim(d)
  if (missing(p)) p <- n
  eig <- purrr::pmap(list(n = n, q = q, p = p), eigen_dec)
  U_eig <- purrr::map(eig, "U")
  U <- U_eig |> rev() |> purrr::reduce(kronecker)
  s_eig <- purrr::map(eig, "s")
  s <- list(rep(s_eig[[1]], p[[2]]), rep(s_eig[[2]], each = p[[1]]))

  sum_wt <- if (reg) sum(wt) else sum(d)
  which_pos <- if (reg) which(wt != 0) else which(ec != 0)
  n_pos <- length(which_pos)
  U_pos <- U[which_pos,]

  if (reg) {

    z <- y
    wt_pos <- c(wt)[which_pos]
    z_pos <- c(z)[which_pos]
    tUWU <- t(U_pos) %*% (wt_pos * U_pos)
    tUWz <- t(U_pos) %*% (wt_pos * z_pos)

  } else {

    off <- log(pmax(ec, 1e-4))
    y <- ifelse(d == 0, NA, log(d)) - off
    y_hat <- log(pmax(d, 1e-8)) - off
    new_wt <- exp(y_hat + off)
  }

  dev_pen <- Inf
  cond_dev_pen <- TRUE

  # Loop
  while (cond_dev_pen) {

    old_dev_pen <- dev_pen

    # update of parameters, working vector and weight matrix
    if (!reg) {

      wt <- new_wt
      z <- y_hat + d / wt - 1
      wt_pos <- c(wt)[which_pos]
      z_pos <- c(z)[which_pos]
      tUWU <- t(U_pos) %*% (wt_pos * U_pos)
      tUWz <- t(U_pos) %*% (wt_pos * z_pos)
    }

    WH_2d_aux <- function(log_lambda) {

      lambda <- exp(log_lambda)

      s_lambda <- purrr::map2(lambda, s, `*`)
      sum_s_lambda <- s_lambda |> do.call(what = `+`)

      Psi_chol <- tUWU
      diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
      Psi_chol <- Psi_chol |> chol()
      Psi <- Psi_chol |> chol2inv()

      beta_hat <- c(Psi %*% tUWz) # fitted value
      y_hat <- c(U %*% beta_hat)
      res <- sqrt(wt) * ((if (reg) y else z) - y_hat) # (weighted) residuals
      dev <- sum(res * res)

      if (criterion != "REML") sum_edf <- sum(Psi * tUWU) # effective degrees of freedom

      score <- switch(criterion,
                      AIC = dev + 2 * sum_edf,
                      BIC = dev + log(prod(n_pos)) * sum_edf,
                      GCV = prod(n_pos) * dev / (prod(n_pos) - sum_edf) ^ 2,
                      REML = {
                        tr_log_P <- purrr::map2(lambda, s, `*`) |> do.call(what = `+`) |> Filter(f = \(x) x > 0) |> log() |> sum()
                        tr_log_Psi <- 2 * (Psi_chol |> diag() |> log() |> sum())
                        RESS <- purrr::map_dbl(s, \(x) sum(beta_hat * x * beta_hat))
                        pen <- purrr::map2(lambda, RESS, `*`) |> do.call(what = `+`)
                        REML <- dev + pen - tr_log_P + tr_log_Psi
                      })

      if (verbose) cat(criterion, " :", format(
        if (criterion == "REML") - score / 2 else score, digits = 3, decimal.mark = ","), "\n")

      return(score)
    }

    lambda <- exp(stats::optim(par = log(lambda),
                               fn = WH_2d_aux,
                               control = list(reltol = accu_crit * sum_wt))$par)

    s_lambda <- purrr::map2(lambda, s, `*`)
    sum_s_lambda <- s_lambda |> do.call(what = `+`)

    Psi_chol <- tUWU
    diag(Psi_chol) <- diag(Psi_chol) + sum_s_lambda
    Psi_chol <- Psi_chol |> chol()
    Psi <- Psi_chol |> chol2inv()

    beta_hat <- c(Psi %*% tUWz) # fitted value

    RESS <- purrr::map_dbl(s, \(x) sum(beta_hat * x * beta_hat))
    edf_par <- colSums(Psi * tUWU) # effective degrees of freedom by parameter
    omega_j <- purrr::map(s_lambda, \(x) ifelse(x == 0, 0, x / sum_s_lambda))

    y_hat <- c(U %*% beta_hat)
    if (!reg) new_wt <- exp(y_hat + off)

    # update of convergence check
    old_dev_pen <- dev_pen

    res <- if (reg) sqrt(wt) * (y - y_hat) else compute_res_deviance(d, new_wt) # (weighted) residuals
    dev <- sum(res * res)
    pen <- purrr::map2(lambda, RESS, `*`) |> do.call(what = `+`)
    dev_pen <- dev + pen

    if (verbose) cat("dev_pen :", format(old_dev_pen, digits = 3, decimal.mark = ","),
                     "=>", format(dev_pen, digits = 3, decimal.mark = ","), "\n")
    cond_dev_pen <- if (reg) FALSE else (old_dev_pen - dev_pen) > accu_dev * sum_wt
  }

  out <- WH_2d_fixed_lambda(d = d, ec = ec, y = y, wt = wt, lambda = lambda, p = n, q = q, reg = reg,
                            verbose = verbose, accu_dev = accu_dev)

  return(out)
}
