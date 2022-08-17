#' WH : A Modern Take on Whittaker-Henderson Smoothing
#'
#' This package implements the popular Whittaker-Henderson Smoothing relying on the latest developments.
#'
#' @import purrr
#'
#' @docType package
#' @name WH-package
NULL

#' Simulated Mortality Datasets
#'
#' Agregated mortality datasets built from simulated portfolios
#'
#' @format A list of 3 agregated datasets obtained from fictive annuity
#'   portfolios of increasing sizes, each one containing :
#'   \describe{\item{exit}{A vector containing the portfolio number of observed
#'   deaths for each age} \item{expo}{A vector containing the portfolio central
#'   exposure in person-years for each age}}
"portfolios_mort"

#' Simulated Long-Term Care Datasets
#'
#' Agregated long-term care datasets built from simulated portfolios
#'
#' @format A list of 3 agregated datasets obtained from fictive LTC annuity
#'   portfolios of increasing sizes, each one containing :
#'   \describe{\item{exit}{A matrix containing the portfolio number of observed
#'   deaths for each combination of age and duration in LTC} \item{expo}{A
#'   matrix containing the portfolio central exposure in person-years for each
#'   combination of age and duration in LTC}}
"portfolios_LTC"

## quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table
# if (getRversion() >= "4.0")  {
#   gv <- c()
#   utils::globalVariables(gv)
# }

# library(devtools)

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
