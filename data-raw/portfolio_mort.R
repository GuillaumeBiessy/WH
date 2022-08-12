library(tidyverse)
library(lubridate)
library(BEAST)

date_tronc_gauche <- dmy("01/01/2008")
date_cens_droite <- dmy("01/01/2018")
causes_sorties <- c("Décès", "Résiliation", "En cours")

portfolio_mort <- read_csv2("data-raw/portefeuille_deces_1M.csv") %>%
  mutate(Sexe = factor(Sexe, levels = c("Homme", "Femme")),
         Date_fin_obs = pmin(Date_deces, Date_resiliation, date_cens_droite, na.rm = TRUE),
         Cause_fin_obs = case_when(Date_fin_obs == Date_deces ~ "Décès",
                                   Date_fin_obs == Date_resiliation ~ "Résiliation",
                                   TRUE ~ "En cours") %>% factor(levels = causes_sorties)) %>%
  select(Date_naissance, Date_souscription, Date_fin_obs, Cause_fin_obs, Sexe) |>
  read_portfolio(age_origin = "Date_naissance",
                 duration_origin = "Date_souscription",
                 obs_start = "Date_souscription",
                 obs_end = "Date_fin_obs",
                 cause = "Cause_fin_obs")

set.seed(1)
portfolios_mort <- c(2e4, 1e5, 5e5) |>
  map(\(x) portfolio_mort[sample(nrow(portfolio_mort), x),]) |>
set_names(c("20k", "100k", "500k")) |>
  map(slice_portfolio,
      timescales = list(age = list(y = 0:120)),
      covariates = list(Sexe = NULL),
      exits = 1) |>
  map(\(x) x[c("expo", "exit")])

usethis::use_data(portfolios_mort, overwrite = TRUE)
