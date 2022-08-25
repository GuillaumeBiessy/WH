data("portfolio_mort")
if(!interactive()) pdf(NULL)

# Data----
keep <- which(portfolio_mort$ec > 0)
d <- portfolio_mort$d[keep]
ec <- portfolio_mort$ec[keep]

y <- log(d / ec)
y[d == 0] <- - 20
wt <- d

# Regression----

## Various way of calling regresion work and method with fixed lambda as well----
ref_fixed_lambda <- WH_1d_reg_fixed_lambda(y, wt, lambda = 1e2)
expect_equal(WH_1d(y = y, wt = wt, lambda = 1e2), ref_fixed_lambda)
expect_equal(WH_1d(d, ec, framework = "reg", lambda = 1e2), ref_fixed_lambda)
expect_equal(WH_1d(d, y = y, lambda = 1e2), ref_fixed_lambda)

## FS is the default method and call FS----
ref_fs <- WH_1d_reg_fs(y, wt)
expect_equal(WH_1d(y = y, wt = wt, method = "fs"), ref_fs)
expect_equal(WH_1d(y = y, wt = wt), ref_fs)

## optim method call optim----
ref_optim <- WH_1d_reg_optim(y, wt)
expect_equal(WH_1d(y = y, wt = wt, method = "optim"), ref_optim)

## optim and fs method match for regression----
expect_equal(ref_fs, ref_optim, tolerance = 1e-6)

## REML is default criterion for optim----
expect_equal(WH_1d(y = y, wt = wt, method = "optim", criterion = "REML"), ref_optim)

## other criteria work----
expect_equal(WH_1d(y = y, wt = wt, criterion = "AIC"),
             WH_1d_reg_optim(y, wt, criterion = "AIC"))
expect_equal(WH_1d(y = y, wt = wt, criterion = "BIC"),
             WH_1d_reg_optim(y, wt, criterion = "BIC"))
expect_equal(WH_1d(y = y, wt = wt, criterion = "GCV"),
             WH_1d_reg_optim(y, wt, criterion = "GCV"))

## rank reduction works----
ref_fs_red <- WH_1d_reg_fs(y, wt, p = 20)
ref_optim_red <- WH_1d_reg_optim(y, wt, p = 20)

expect_equal(WH_1d(y = y, wt = wt, p = 20), ref_fs_red)
expect_equal(WH_1d(y = y, wt = wt, method = "optim", p = 20), ref_optim_red)
expect_equal(ref_fs_red, ref_optim_red, tolerance = 1e-6)

# Maximum likelihood----

## fixed lambda method works----
expect_equal(WH_1d(d, ec, lambda = 1e2),
             WH_1d_ml_fixed_lambda(d, ec, lambda = 1e2))

## FS method works and is default----
ref_ml_fs <- WH_1d_ml_fs(d, ec)
expect_equal(WH_1d(d, ec), ref_ml_fs)

## optim method works----
ref_ml_optim <- WH_1d_ml_optim(d, ec)
expect_equal(WH_1d(d, ec, method = "optim"),
             ref_ml_optim)
expect_equal(WH_1d_ml_optim(d, ec, criterion = "REML"),
             ref_ml_optim)

## optim and fs method are not too far away for ML----
expect_equal(ref_ml_fs, ref_ml_optim, tolerance = 1e-1)

## other criteria work----
expect_equal(WH_1d(d, ec, criterion = "AIC"),
             WH_1d_ml_optim(d, ec, criterion = "AIC"))
expect_equal(WH_1d(d, ec, criterion = "BIC"),
             WH_1d_ml_optim(d, ec, criterion = "BIC"))
expect_equal(WH_1d(d, ec, criterion = "GCV"),
             WH_1d_ml_optim(d, ec, criterion = "GCV"))

## rank reduction works----
ref_fs_red <- WH_1d_ml_fs(d, ec, p = 20)
ref_optim_red <- WH_1d_ml_optim(d, ec, p = 20)
expect_equal(WH_1d(d, ec, p = 20), ref_fs_red)
expect_equal(WH_1d(d, ec, method = "optim", p = 20), ref_optim_red)
expect_equal(ref_fs_red, ref_optim_red, tolerance = 1e-1)

# Extrapolation----

newdata <- 18:99

expect_equal(WH_1d_reg_fs(y, wt) |> predict(newdata),
             WH_1d_reg_optim(y, wt) |> predict(newdata),
             tolerance = 1e-6)

expect_equal(WH_1d_ml_fs(d, ec) |> predict(newdata),
             WH_1d_ml_optim(d, ec) |> predict(newdata),
             tolerance = 1e-1)

## Regression----
fs_extra_reg <- WH_1d_reg_fs(y, wt) |> predict(newdata)
optim_extra_reg <- WH_1d_reg_fs(y, wt) |> predict(newdata)
expect_equal(fs_extra_reg, optim_extra_reg, tolerance = 1e-6)

# Maximum likelihood----
fs_extra_ml <- WH_1d_reg_fs(y, wt) |> predict(newdata)
optim_extra_ml <- WH_1d_reg_fs(y, wt) |> predict(newdata)
expect_equal(fs_extra_ml, optim_extra_ml, tolerance = 1e-1)

# Plots----

## Regression----
WH_1d_reg_optim(y, wt) |> plot()
WH_1d_reg_fs(y, wt) |> plot()

WH_1d_reg_optim(y, wt) |> plot("res")
WH_1d_reg_fs(y, wt) |> plot("res")

WH_1d_reg_optim(y, wt) |> plot("edf")
WH_1d_reg_fs(y, wt) |> plot("edf")

# Maximum likelihood----
WH_1d_ml_optim(d, ec) |> plot()
WH_1d_ml_fs(d, ec) |> plot()

WH_1d_ml_optim(d, ec) |> plot("res")
WH_1d_ml_fs(d, ec) |> plot("res")

WH_1d_ml_optim(d, ec) |> plot("edf")
WH_1d_ml_fs(d, ec) |> plot("edf")

## Extrapolation----
fs_extra_reg |> plot()
optim_extra_reg |> plot()

fs_extra_ml |> plot()
optim_extra_ml |> plot()

