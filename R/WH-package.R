#' WH : A Modern Take on Whittaker-Henderson Smoothing
#'
#' Implementation of the popular Whittaker-Henderson smoothing method for the
#' gradation of actuarial tables for person insurance risks such as mortality,
#' disability, long-term care and unemployment. The package handles tables in 1
#' or 2 dimensions and performs automatic selection of the smoothing parameters.
#' Smoothing is applied to vectors / matrices containing the number of observed
#' events and the associated central risk exposure in the framework of maximum
#' likelihood estimation or an approximate regression framework. Optimal
#' parameter selection relies on [stats::optimize()] function and the
#' Nelder-Mead algorithm from the [stats::optim()] function in the stats package
#' or uses the generalized Fellner-Schall method.
#'
#' @docType package
#' @name WH-package
NULL

#' Agregated Mortality Dataset
#'
#' Agregated dataset built from a simulated mortality portfolio
#'
#' @format An agregated dataset agregating the information from an annuity
#'   portfolio with 100,000 contributors on a 10-year observation period. The
#'   dataset is supplied as a list with two components : \describe{\item{exit}{A
#'   vector containing the portfolio number of observed deaths for each age from
#'   0 to 120} \item{expo}{A vector containing the portfolio central exposure in
#'   person-years for each age from 0 to 119}}
"portfolio_mort"

#' Agregated Long-Term Care Dataset
#'
#' Agregated dataset built from a simulated long-term care portfolio
#'
#' @format An agregated dataset obtaining by agregating the information from a
#'   fictive annuitant database from a long-term care cover with 5,000
#'   annuitants on a 10-year observation period. The dataset is supplied as a
#'   list with two components : \describe{\item{exit}{A matrix
#'   containing the portfolio number of observed deaths for each combination of
#'   age from 0 to 119 and duration in LTC from 0 to 29} \item{expo}{A matrix
#'   containing the portfolio central exposure in person-years for each
#'   combination of age from 0 to 119 and duration in LTC from 0 to 29}}
"portfolio_LTC"

# Import----
# use_package("purrr")

# Data----
# use_data_raw("portfolio_mort")
# use_data_raw("portfolio_LTC")

# Tests----
# use_test("1D")
# use_test("2D")

# Patchnotes----
# use_news_md()
