library(tidyverse)
library(lubridate)
library(BEAST)

date_tronc_gauche <- dmy("01/01/2008")
date_cens_droite <- dmy("01/01/2018")
franchise <- days(90)
causes_sorties <- c("Décès", "En cours")

portfolio_LTC <- read_csv2("data-raw/portefeuille_dependance_1M.csv") |>
  filter(!is.na(Date_dependance)) |>
  mutate(Sexe = factor(Sexe, levels = c("Homme", "Femme")),
         Date_debut_obs = pmax(Date_dependance + franchise, date_tronc_gauche, na.rm = TRUE),
         Date_fin_obs = pmin(Date_deces, date_cens_droite, na.rm = TRUE),
         Cause_fin_obs = case_when(Date_fin_obs == Date_deces ~ "Décès",
                                   TRUE ~ "En cours") %>% factor(levels = causes_sorties)) |>
  select(Date_naissance, Date_dependance, Date_debut_obs, Date_fin_obs, Cause_fin_obs, Sexe) |>
  read_portfolio(age_origin = "Date_naissance",
                 duration_origin = "Date_dependance",
                 obs_start = "Date_debut_obs",
                 obs_end = "Date_fin_obs",
                 cause = "Cause_fin_obs")

set.seed(1)
portfolio_LTC <- c(1e3, 5e3, 2.5e4) |>
  map(\(x) portfolio_LTC[sample(nrow(portfolio_LTC), x),]) |>
  set_names(c("1k", "5k", "25k")) |>
  map(slice_portfolio,
      timescales = list(age = list(y = 0:120),
                        duration = list(y = seq(0, 30, 1))),
      covariates = list(Sexe = NULL),
      exits = 1) |>
  map(\(x) list(d = x$exit |> aperm(c(3,4,1,2)) |> colSums(dims = 2),
                ec = (x$expo / 365.25)  |> aperm(c(3,1,2)) |> colSums())) |>
  (\(x) x[[1]])()

usethis::use_data(portfolio_LTC, overwrite = TRUE)

# set.seed(1)
# portfolios_LTC <- c(1e3, 5e3, 2.5e4) |>
#   map(\(x) portfolio_LTC[sample(nrow(portfolio_LTC), x),]) |>
#   set_names(c("1k", "5k", "25k")) |>
#   map(slice_portfolio,
#       timescales = list(age = list(y = 0:120),
#                         duration = list(y = seq(0, 30, 1))),
#       covariates = list(Sexe = NULL),
#       exits = 1) |>
#   map(\(x) list(d = x$exit |> aperm(c(3,4,1,2)) |> colSums(dims = 2),
#                 ec = (x$expo / 365.25)  |> aperm(c(3,1,2)) |> colSums()))
#
# usethis::use_data(portfolios_LTC, overwrite = TRUE)
