library(tidyverse)
library(lubridate)

date_tronc_gauche <- dmy("01/01/2008")
date_cens_droite <- dmy("01/01/2018")
causes_sorties <- c("Décès", "Résiliation", "En cours")

portfolio_mort <- read_csv2("data-raw/portefeuille_deces.csv") %>%
  mutate(Sexe = factor(Sexe, levels = c("Homme", "Femme")),
         Date_fin_obs = pmin(Date_deces, Date_resiliation, date_cens_droite, na.rm = TRUE),
         Cause_fin_obs = case_when(Date_fin_obs == Date_deces ~ "Décès",
                                   Date_fin_obs == Date_resiliation ~ "Résiliation",
                                   TRUE ~ "En cours") %>% factor(levels = causes_sorties)) %>%
  select(Date_naissance, Date_souscription, Date_fin_obs, Cause_fin_obs, Sexe)

usethis::use_data(portfolio_mort, overwrite = TRUE)
