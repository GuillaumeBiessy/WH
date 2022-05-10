data("portfolio_mort")
data("portfolio_LTC")

library(profvis)
library(bench)

# 1D smoothing

type_IC <- "freq"

ec <- rowSums(portfolio_mort$expo) / 365.25
d <- rowSums(portfolio_mort$exit)

y <- log(d / ec)
y[d == 0] <- 0

WH_1d_fit_reg_vintage <- WH_1d_reg_vintage(y, d, type_IC = type_IC)
WH_1d_fit_reg_mixed <- WH_1d_reg_mixed(y, d, type_IC = type_IC)

WH_1d_fit_surv_vintage <- WH_1d_surv_vintage(ec, d, type_IC = type_IC)
WH_1d_fit_surv_mixed <- WH_1d_surv_mixed(ec, d, type_IC = type_IC)

plot.WH_1d(WH_1d_fit_reg_vintage)
plot.WH_1d(WH_1d_fit_reg_mixed)
plot.WH_1d(WH_1d_fit_surv_vintage)
plot.WH_1d(WH_1d_fit_surv_mixed)

newdata_1d <- 20:109

WH_1d_fit_reg_vintage_pred <- predict.WH_1d(object = WH_1d_fit_reg_vintage,
                                            newdata = newdata_1d, type_IC = type_IC)
WH_1d_fit_reg_mixed_pred <- predict.WH_1d(object = WH_1d_fit_reg_mixed,
                                          newdata = newdata_1d, type_IC = type_IC)
WH_1d_fit_surv_vintage_pred <- predict.WH_1d(object = WH_1d_fit_surv_vintage,
                                             newdata = newdata_1d, type_IC = type_IC)
WH_1d_fit_surv_mixed_pred <- predict.WH_1d(object = WH_1d_fit_surv_mixed,
                                           newdata = newdata_1d, type_IC = type_IC)

plot.WH_1d(WH_1d_fit_reg_vintage_pred)
plot.WH_1d(WH_1d_fit_reg_mixed_pred)
plot.WH_1d(WH_1d_fit_surv_vintage_pred)
plot.WH_1d(WH_1d_fit_surv_mixed_pred)

library(ggplot2)
library(Linkplots)

theme_LinkPact() |> theme_set()

# 2D smoothing

EC <- (portfolio_LTC$expo / 365.25) |> aperm(c(3,1,2)) |> colSums()
D <- (portfolio_LTC$exit) |> aperm(c(3,4,1,2)) |> colSums(dims = 2)

Y <- log(D / EC) # observation vector
Y[D == 0] <- 0

WH_2d_fit_reg_vintage <- WH_2d_reg_vintage(Y, D, type_IC = type_IC)
WH_2d_fit_reg_mixed <- WH_2d_reg_mixed(Y, D, type_IC = type_IC)

WH_2d_fit_surv_vintage <- WH_2d_surv_vintage(EC, D, type_IC = type_IC)
WH_2d_fit_surv_mixed <- WH_2d_surv_mixed(EC, D, type_IC = type_IC)

plot.WH_2d(WH_2d_fit_reg_vintage)
plot.WH_2d(WH_2d_fit_reg_mixed)
plot.WH_2d(WH_2d_fit_surv_vintage)
plot.WH_2d(WH_2d_fit_surv_mixed)

newdata_2d <- list(age = 55:109,
                duration = 0:19)

WH_2d_fit_reg_vintage_pred <- predict.WH_2d(object = WH_2d_fit_reg_vintage,
                                            newdata = newdata_2d, type_IC = type_IC)
WH_2d_fit_reg_mixed_pred <- predict.WH_2d(object = WH_2d_fit_reg_mixed,
                                          newdata = newdata_2d, type_IC = type_IC)
WH_2d_fit_surv_vintage_pred <- predict.WH_2d(object = WH_2d_fit_surv_vintage,
                                             newdata = newdata_2d, type_IC = type_IC)
WH_2d_fit_surv_mixed_pred <- predict.WH_2d(object = WH_2d_fit_surv_mixed,
                                           newdata = newdata_2d, type_IC = type_IC)

plot.WH_2d(WH_2d_fit_reg_vintage_pred)
plot.WH_2d(WH_2d_fit_reg_mixed_pred)
plot.WH_2d(WH_2d_fit_surv_vintage_pred)
plot.WH_2d(WH_2d_fit_surv_mixed_pred)
