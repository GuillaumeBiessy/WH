% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{predict.WH_2d}
\alias{predict.WH_2d}
\title{Prediction for a Whittaker-Henderson Fit}
\usage{
\method{predict}{WH_2d}(object, newdata = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{"WH_2d"} returned by the \code{\link[=WH_2d]{WH_2d()}} function}

\item{newdata}{A list containing two vectors indicating the new observation
positions}

\item{...}{Not used}
}
\value{
An object of class \code{"WH_2d"} with additional components \code{y_pred} and
\code{std_y_pred} corresponding to the model predictions and associated standard
deviations.
}
\description{
Extrapolate the Whittaker-Henderson fit for new observations in a way that is
consistent with the initial model fit.
}
\examples{
keep_age <- which(rowSums(portfolio_LTC$ec) > 5e2)
keep_duration <- which(colSums(portfolio_LTC$ec) > 1e3)

d  <- portfolio_LTC$d[keep_age, keep_duration]
ec <- portfolio_LTC$ec[keep_age, keep_duration]

fit <- WH_2d(d, ec)
newdata <- list(age = 50:99, duration = 0:19)
pred <- predict(fit, newdata)
plot(pred)

}
