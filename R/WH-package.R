#' WH : Enhanced Implementation of Whittaker-Henderson Smoothing
#'
#' An enhanced implementation of Whittaker-Henderson smoothing for the gradation
#' of one-dimensional and two-dimensional actuarial tables used to quantify Life
#' Insurance risks. `WH` is based on the methods described in Biessy (2023)
#' <arXiv:2306.06932>. Among other features, it generalizes the original
#' smoothing algorithm to maximum likelihood estimation, automatically selects
#' the smoothing parameter(s) and extrapolates beyond the range of data.
#'
#' @docType package
#' @name WH-package
NULL

#' Agregated Mortality Dataset
#'
#' Agregated dataset built from a synthetic mortality portfolio
#'
#' @format A dataset containing the information from a simulated annuity
#'   portfolio with 100,000 contributors over a 10-year observation period. The
#'   dataset is supplied as a list with two components : \describe{\item{d}{A
#'   vector containing the number of observed deaths for ages where at
#'   least one death has been observed} \item{ec}{A vector containing the associated
#'   central exposure in person-years for each age in d}}
"portfolio_mort"

#' Agregated Long-Term Care Dataset
#'
#' Agregated dataset built from a synthetic long-term care portfolio
#'
#' @format A dataset obtaining by agregating the information from a fictive
#'   long-term care annuitant database with 5,000 annuitants over a 10-year
#'   observation period. The dataset is supplied as a list with two components :
#'   \describe{\item{exit}{A matrix containing the number of observed deaths for
#'   each combination of age and duration in LTC where at least one death
#'   has been observed} \item{expo}{A matrix containing the associated central
#'   exposure in person-years for each combination of age and duration in LTC in d}}
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
