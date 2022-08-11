data("portfolios_mort")
data("portfolios_LTC")

get_edf <- \(x) sum(x$edf)

# 1D smoothing----
ec <- rowSums(portfolios_mort[[1]]$expo) / 365.25
d <- rowSums(portfolios_mort[[1]]$exit)

keep <- which(ec > 0) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()
ec <- ec[keep]
d <- d[keep]

y <- log(d / ec)
y[d == 0] <- - 20
wt <- d

# Regression
WH_1d_reg_fixed_lambda(y, wt, lambda = 1e2) |> get_edf()
WH_1d(y = y, wt = wt, lambda = 1e2) |> get_edf()
WH_1d(d, ec, framework = "reg", lambda = 1e2) |> get_edf()

WH_1d_reg_fs(y, wt) |> get_edf()
WH_1d(y = y, wt = wt) |> get_edf()
WH_1d(d, ec, framework = "reg") |> get_edf()

WH_1d_reg_optim(y, wt) |> get_edf()
WH_1d(y = y, wt = wt, method = "optim") |> get_edf()
WH_1d(d, ec, framework = "reg", method = "optim") |> get_edf()

WH_1d_reg_optim(y, wt, criterion = "AIC") |> get_edf()
WH_1d(y = y, wt = wt, method = "optim", criterion = "AIC") |> get_edf()
WH_1d(d, ec, framework = "reg", method = "optim", criterion = "AIC") |> get_edf()

WH_1d_reg_optim(y, wt, criterion = "BIC") |> get_edf()
WH_1d(y = y, wt = wt, method = "optim", criterion = "BIC") |> get_edf()
WH_1d(d, ec, framework = "reg", method = "optim", criterion = "BIC") |> get_edf()

WH_1d_reg_optim(y, wt, criterion = "GCV") |> get_edf()
WH_1d(y = y, wt = wt, method = "optim", criterion = "GCV") |> get_edf()
WH_1d(d, ec, framework = "reg", method = "optim", criterion = "GCV") |> get_edf()

WH_1d_reg_fs(y, wt, p = 20) |> get_edf()
WH_1d(y = y, wt = wt, p = 20) |> get_edf()
WH_1d(d, ec, framework = "reg", p = 20) |> get_edf()

WH_1d_reg_optim(y, wt, p = 20) |> get_edf()
WH_1d(y = y, wt = wt, method = "optim", p = 20) |> get_edf()
WH_1d(d, ec, framework = "reg", method = "optim", p = 20) |> get_edf()

# Maximum likelihood
WH_1d_ml_fixed_lambda(d, ec, lambda = 1e2) |> get_edf()

WH_1d_ml_fs(d, ec) |> get_edf()

WH_1d_ml_optim(d, ec) |> get_edf()

WH_1d_ml_optim(d, ec, criterion = "AIC") |> get_edf()
WH_1d_ml_optim(d, ec, criterion = "BIC") |> get_edf()
WH_1d_ml_optim(d, ec, criterion = "GCV") |> get_edf()

WH_1d(d, ec, p = 20) |> get_edf()

WH_1d(d, ec, method = "optim", p = 20) |> get_edf()

# Plots
WH_1d_reg_optim(y, wt) |> plot.WH_1d()
WH_1d_reg_fs(y, wt) |> plot.WH_1d()

WH_1d_ml_optim(d, ec) |> plot.WH_1d()
WH_1d_ml_fs(d, ec) |> plot.WH_1d()

WH_1d_reg_optim(y, wt) |> plot.WH_1d("res")
WH_1d_reg_fs(y, wt) |> plot.WH_1d("res")

WH_1d_ml_optim(d, ec) |> plot.WH_1d("res")
WH_1d_ml_fs(d, ec) |> plot.WH_1d("res")

WH_1d_reg_optim(y, wt) |> plot.WH_1d("edf")
WH_1d_reg_fs(y, wt) |> plot.WH_1d("edf")

WH_1d_ml_optim(d, ec) |> plot.WH_1d("edf")
WH_1d_ml_fs(d, ec) |> plot.WH_1d("edf")

# Extrapolation
newdata_1d <- 18:119

WH_1d_reg_optim(y, wt) |> predict.WH_1d(newdata_1d) |> plot.WH_1d()
WH_1d_reg_fs(y, wt) |> predict.WH_1d(newdata_1d) |> plot.WH_1d()

WH_1d_ml_optim(d, ec) |> predict.WH_1d(newdata_1d) |> plot.WH_1d()
WH_1d_ml_fs(d, ec) |> predict.WH_1d(newdata_1d) |> plot.WH_1d()

# 2D smoothing----
d  <- (portfolios_LTC[[1]]$exit) |> aperm(c(3,4,1,2)) |> colSums(dims = 2)
ec <- (portfolios_LTC[[1]]$expo / 365.25) |> aperm(c(3,1,2)) |> colSums()

keep_age <- which(rowSums(ec) > 1e2) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()
keep_duration <- which(colSums(ec) > 1e2) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()

d <- d[keep_age, keep_duration]
ec <- ec[keep_age, keep_duration]

y <- log(d / ec) # observation vector
y[d == 0] <- - 20
wt <- d

# Regression
WH_2d_reg_fixed_lambda(y, wt, lambda = c(1e2, 1e2)) |> get_edf()
WH_2d(y = y, wt = wt, lambda = c(1e2, 1e2)) |> get_edf()
WH_2d(d, ec, framework = "reg", lambda = c(1e2, 1e2)) |> get_edf()

WH_2d_reg_fs(y, wt) |> get_edf()
WH_2d(y = y, wt = wt) |> get_edf()
WH_2d(d, ec, framework = "reg") |> get_edf()

WH_2d_reg_optim(y, wt) |> get_edf()
WH_2d(y = y, wt = wt, method = "optim") |> get_edf()
WH_2d(d, ec, framework = "reg", method = "optim") |> get_edf()

WH_2d_reg_optim(y, wt, criterion = "AIC") |> get_edf()
WH_2d(y = y, wt = wt, method = "optim", criterion = "AIC") |> get_edf()
WH_2d(d, ec, framework = "reg", method = "optim", criterion = "AIC") |> get_edf()

WH_2d_reg_optim(y, wt, criterion = "BIC") |> get_edf()
WH_2d(y = y, wt = wt, method = "optim", criterion = "BIC") |> get_edf()
WH_2d(d, ec, framework = "reg", method = "optim", criterion = "BIC") |> get_edf()

WH_2d_reg_optim(y, wt, criterion = "GCV") |> get_edf()
WH_2d(y = y, wt = wt, method = "optim", criterion = "GCV") |> get_edf()
WH_2d(d, ec, framework = "reg", method = "optim", criterion = "GCV") |> get_edf()

WH_2d_reg_fs(y, wt, p = c(10, 5)) |> get_edf()
WH_2d(y = y, wt = wt, p = c(10, 5)) |> get_edf()
WH_2d(d, ec, framework = "reg", p = c(10, 5)) |> get_edf()

WH_2d_reg_optim(y, wt, p = c(10, 5)) |> get_edf()
WH_2d(y = y, wt = wt, method = "optim", p = c(10, 5)) |> get_edf()
WH_2d(d, ec, framework = "reg", method = "optim", p = c(10, 5)) |> get_edf()

# Maximum likelihood
WH_2d_ml_fixed_lambda(d, ec, lambda = c(1e2, 1e2)) |> get_edf()

WH_2d_ml_fs(d, ec) |> get_edf()

WH_2d_ml_optim(d, ec) |> get_edf()
WH_2d_ml_optim(d, ec, criterion = "AIC") |> get_edf()
WH_2d_ml_optim(d, ec, criterion = "BIC") |> get_edf()
WH_2d_ml_optim(d, ec, criterion = "GCV") |> get_edf()

WH_2d(d, ec, p = c(10, 5)) |> get_edf()
WH_2d(d, ec, method = "optim", p = c(10, 5)) |> get_edf()

# Plots
WH_2d_reg_optim(y, wt) |> plot.WH_2d()
WH_2d_reg_fs(y, wt) |> plot.WH_2d()

WH_2d_ml_optim(d, ec) |> plot.WH_2d()
WH_2d_ml_fs(d, ec) |> plot.WH_2d()

WH_2d_reg_optim(y, wt) |> plot.WH_2d("std_y_hat")
WH_2d_reg_fs(y, wt) |> plot.WH_2d("std_y_hat")

WH_2d_ml_optim(d, ec) |> plot.WH_2d("std_y_hat")
WH_2d_ml_fs(d, ec) |> plot.WH_2d("std_y_hat")

# Extrapolation
newdata_2d <- list(age = 50:119,
                duration = 0:29)

WH_2d_reg_optim(y, wt) |> predict.WH_2d(newdata_2d) |> plot.WH_2d()
WH_2d_reg_fs(y, wt) |> predict.WH_2d(newdata_2d) |> plot.WH_2d()

WH_2d_ml_optim(d, ec) |> predict.WH_2d(newdata_2d) |> plot.WH_2d()
WH_2d_ml_fs(d, ec) |> predict.WH_2d(newdata_2d) |> plot.WH_2d()

