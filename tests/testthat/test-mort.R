data("portfolio_mort")
data("portfolio_LTC")

library(BEAST)
library(profvis)
library(bench)
library(arrays)

portfolio_mort_agrege <- portfolio_mort |>
  read_portfolio(age_origin = "Date_naissance",
                 duration_origin = "Date_souscription",
                 obs_start = "Date_souscription",
                 obs_end = "Date_fin_obs",
                 cause = "Cause_fin_obs") |>
  slice_portfolio(timescales = list(age = list(y = 30:90)),
                  covariates = list(Sexe = NULL),
                  exits = 1)

portfolio_LTC_agrege <- portfolio_LTC |>
  read_portfolio(age_origin = "Date_naissance",
                 duration_origin = "Date_dependance",
                 obs_start = "Date_debut_obs",
                 obs_end = "Date_fin_obs",
                 cause = "Cause_fin_obs") |>
  slice_portfolio(timescales = list(age = list(y = 65:100),
                                    duration = list(m = seq(1, 120, 3))),
                  covariates = list(Sexe = NULL),
                  exits = 1)

E <- rowSums(portfolio_mort_agrege$expo) / 365.25
D <- rowSums(portfolio_mort_agrege$exit)

WH_1d_vintage(y, wt)
WH_1d(y, wt) |> plot.WH_1d()

E <- aSums(portfolio_LTC_agrege$expo, 1:2) / 365.25
D <- aSums(portfolio_LTC_agrege$exit, 1:2)

y <- log(D / E) # observation vector
y[D == 0] <- 0
wt <- D

WH_2d_vintage(y, wt)
WH_2d(y, wt) |> plot.WH_2d()
