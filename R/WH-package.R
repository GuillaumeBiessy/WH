#' WH : Enhanced Implementation of Whittaker-Henderson Smoothing
#'
#' An enhanced implementation of Whittaker-Henderson smoothing for the gradation
#' of one-dimensional and two-dimensional actuarial tables used to quantify Life
#' Insurance risks. `WH` is based on the methods described in Biessy (2025)
#' <doi:10.48550/arXiv.2306.06932>. Among other features, it generalizes the
#' original smoothing algorithm to maximum likelihood estimation, automatically
#' selects the smoothing parameter(s) and extrapolates beyond the range of data.
#'
#' @name WH-package
"_PACKAGE"

#' @useDynLib WH, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Agregated Mortality Dataset
#'
#' Agregated dataset built from a simulated mortality portfolio
#'
#' @format An synthetic aggregated dataset with death and exposure counts from a
#'   simulated annuity portfolio with 100,000 contributors on a 20-year
#'   observation period. The
#'   dataset is supplied as a list with two components : \describe{\item{d}{A
#'   vector containing the portfolio number of observed deaths for each age from
#'   50 to 95 (excluded)} \item{ec}{A vector containing the portfolio central exposure in
#'   person-years for each age from 50 to 95 (excluded)}}
"portfolio_mort"

#' Agregated Long-Term Care Dataset
#'
#' Agregated dataset built from a simulated long-term care portfolio
#'
#' @format An synthetic aggregated dataset with death and exposure counts from a
#'   simulated long-term care portfolio with 100,000 contributors on a 20-year
#'   observation period (only deaths following long-term care claims are
#'   counted). The dataset is supplied as a
#'   list with two components : \describe{\item{d}{A matrix
#'   containing the portfolio number of observed deaths for each combination of
#'   age from 70 to 100 (excluded) and duration in LTC from 0 to 15 (excluded)} \item{ec}{A matrix
#'   containing the portfolio central exposure in person-years for each
#'   combination age from 70 to 100 (excluded) and duration in LTC from 0 to 15 (excluded)}}
"portfolio_LTC"

# Import----
# use_package("stats")

# Data----
# use_data_raw("portfolio_mort")
# use_data_raw("portfolio_LTC")

# Tests----
# use_test("1D")
# use_test("2D")

# Patchnotes----
# use_news_md()
