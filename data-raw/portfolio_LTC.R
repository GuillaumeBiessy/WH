library(tidyverse)
library(lubridate)

date_tronc_gauche <- dmy("01/01/2008")
date_cens_droite <- dmy("01/01/2018")
franchise <- days(90)
causes_sorties <- c("Décès", "En cours")

portfolio_LTC <- read_csv2("data-raw/portefeuille_dependance.csv") %>%
  filter(!is.na(Date_dependance)) %>%
  mutate(Sexe = factor(Sexe, levels = c("Homme", "Femme")),
         Date_debut_obs = pmax(Date_dependance + franchise, date_tronc_gauche, na.rm = TRUE),
         Date_fin_obs = pmin(Date_deces, date_cens_droite, na.rm = TRUE),
         Cause_fin_obs = case_when(Date_fin_obs == Date_deces ~ "Décès",
                                   TRUE ~ "En cours") %>% factor(levels = causes_sorties)) %>%
  select(Date_naissance, Date_dependance, Date_debut_obs, Date_fin_obs, Cause_fin_obs, Sexe)

usethis::use_data(portfolio_LTC, overwrite = TRUE)
