#' WH : a modern take on Whittaker-Henderson smoothing
#'
#' The package implements the popular (at least among actuaries)
#' Whittaker-Henderson smoothing relying on the latest developments.
#'
#' @import purrr
#'
#' @docType package
#' @name WH-package
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table
if (getRversion() >= "4.0")  {
  gv <- c()
  utils::globalVariables(gv)
}

# library(devtools)

# Import----
# use_package("purrr")

# Data----
# use_data_raw("portfolio_mort")
# use_data_raw("portfolio_LTC")

# Tests----
# use_test("mort")

# Patchnotes----
# use_news_md()
