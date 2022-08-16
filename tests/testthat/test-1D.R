data("portfolios_mort")

get_edf <- \(x) sum(x$edf)

# 1D smoothing----
d <- portfolios_mort[[1]]$d
ec <- portfolios_mort[[1]]$ec

keep <- which(ec > 0) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()
d <- d[keep]
ec <- ec[keep]

y <- log(d / ec)
y[d == 0] <- - 20
wt <- d

# Regression

ref_fixed_lambda <- WH_1d_reg_fixed_lambda(y, wt, lambda = 1e2)
expect_equal(WH_1d(y = y, wt = wt, lambda = 1e2), ref_fixed_lambda)
expect_equal(WH_1d(d, ec, framework = "reg", lambda = 1e2), ref_fixed_lambda)

ref_fs <- WH_1d_reg_fs(y, wt)
expect_equal(WH_1d(y = y, wt = wt), ref_fs)
expect_equal(WH_1d(d, ec, framework = "reg"), ref_fs)

ref_optim <- WH_1d_reg_optim(y, wt)
expect_equal(WH_1d(y = y, wt = wt, method = "optim"), ref_optim)
expect_equal(WH_1d(d, ec, framework = "reg", method = "optim"), ref_optim)

expect_equal(ref_fs, ref_optim, tolerance = 1e-6)

ref_aic <- WH_1d_reg_optim(y, wt, criterion = "AIC")
expect_equal(WH_1d(y = y, wt = wt, method = "optim", criterion = "AIC"), ref_aic)
expect_equal(WH_1d(d, ec, framework = "reg", method = "optim", criterion = "AIC"), ref_aic)

ref_bic <- WH_1d_reg_optim(y, wt, criterion = "BIC")
expect_equal(WH_1d(y = y, wt = wt, method = "optim", criterion = "BIC"), ref_bic)
expect_equal(WH_1d(d, ec, framework = "reg", method = "optim", criterion = "BIC"), ref_bic)

ref_gcv <- WH_1d_reg_optim(y, wt, criterion = "GCV")
expect_equal(WH_1d(y = y, wt = wt, method = "optim", criterion = "GCV"), ref_gcv)
expect_equal(WH_1d(d, ec, framework = "reg", method = "optim", criterion = "GCV"), ref_gcv)

ref_fs_red <- WH_1d_reg_fs(y, wt, p = 20)
expect_equal(WH_1d(y = y, wt = wt, p = 20), ref_fs_red)
expect_equal(WH_1d(d, ec, framework = "reg", p = 20), ref_fs_red)

ref_optim_red <- WH_1d_reg_optim(y, wt, p = 20)
expect_equal(WH_1d(y = y, wt = wt, method = "optim", p = 20), ref_optim_red)
expect_equal(WH_1d(d, ec, framework = "reg", method = "optim", p = 20), ref_optim_red)

expect_equal(ref_fs_red, ref_optim_red, tolerance = 1e-6)

# Maximum likelihood
expect_equal(WH_1d(d, ec, lambda = 1e2),
             WH_1d_ml_fixed_lambda(d, ec, lambda = 1e2))

ref_ml_fs <- WH_1d_ml_fs(d, ec)
expect_equal(WH_1d_ml_fs(d, ec), ref_ml_fs)

ref_ml_optim <- WH_1d_ml_optim(d, ec)
expect_equal(WH_1d_ml_optim(d, ec),
             ref_ml_optim)

expect_equal(ref_ml_fs, ref_ml_optim, tolerance = 1e-1)

WH_1d_ml_optim(d, ec, criterion = "AIC")
WH_1d_ml_optim(d, ec, criterion = "BIC")
WH_1d_ml_optim(d, ec, criterion = "GCV")

WH_1d(d, ec, p = 20)

WH_1d(d, ec, method = "optim", p = 20)

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
