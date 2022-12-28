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
ref_fixed_lambda <- WH_1d_fixed_lambda(y = y, wt = wt, lambda = 1e2, reg = TRUE)
expect_equal(WH_1d(y = y, wt = wt, lambda = 1e2), ref_fixed_lambda)
expect_equal(WH_1d(d, ec, framework = "reg", lambda = 1e2), ref_fixed_lambda)
expect_equal(WH_1d(d, y = y, lambda = 1e2), ref_fixed_lambda)

## outer is the default method and calls outer----
ref_outer <- WH_1d_outer(y = y, wt = wt, reg = TRUE)
expect_equal(WH_1d(y = y, wt = wt, method = "outer"), ref_outer)
expect_equal(WH_1d(y = y, wt = wt), ref_outer)

## perf method calls perf----
ref_perf <- WH_1d_perf(y = y, wt = wt, reg = TRUE)
expect_equal(WH_1d(y = y, wt = wt, method = "perf"), ref_perf)

## fs method calls fs----
ref_fs <- WH_1d_fs(y = y, wt = wt, reg = TRUE)
expect_equal(WH_1d(y = y, wt = wt, method = "fs"), ref_fs)

## all methods match for regression----
expect_equal(ref_fs, ref_outer, tolerance = 1e-5)
expect_equal(ref_perf, ref_outer, tolerance = 1e-5)

## outer is default criterion for optim----
expect_equal(WH_1d(y = y, wt = wt, criterion = "REML"), ref_outer)

## other criteria work----
expect_equal(WH_1d(y = y, wt = wt, criterion = "AIC"),
             WH_1d_outer(y = y, wt = wt, criterion = "AIC", reg = TRUE))
expect_equal(WH_1d(y = y, wt = wt, criterion = "BIC"),
             WH_1d_outer(y = y, wt = wt, criterion = "BIC", reg = TRUE))
expect_equal(WH_1d(y = y, wt = wt, criterion = "GCV"),
             WH_1d_outer(y = y, wt = wt, criterion = "GCV", reg = TRUE))

## rank reduction works----
ref_perf_red <- WH_1d_perf(y = y, wt = wt, p = 20, reg = TRUE)
ref_outer_red <- WH_1d_outer(y = y, wt = wt, p = 20, reg = TRUE)

expect_equal(WH_1d(y = y, wt = wt, p = 20), ref_outer_red)
expect_equal(WH_1d(y = y, wt = wt, method = "perf", p = 20), ref_perf_red)
expect_equal(ref_perf_red, ref_outer_red, tolerance = 1e-5)

# Maximum likelihood----

## fixed lambda method works----

expect_equal(WH_1d(d, ec, lambda = 1e2),
             WH_1d_fixed_lambda(d, ec, lambda = 1e2))

## perf method works and is default----
ref_ml_perf <- WH_1d_perf(d, ec)
expect_equal(WH_1d(d, ec, method = "perf"), ref_ml_perf)

## optim method works----
ref_ml_outer <- WH_1d_outer(d, ec)
expect_equal(WH_1d(d, ec, method = "outer"),
             ref_ml_outer)
expect_equal(WH_1d_outer(d, ec, criterion = "REML"),
             ref_ml_outer)

## optim and perf method are not too far away for ML----
expect_equal(ref_ml_perf, ref_ml_outer, tolerance = 1e-1)

## other criteria work----
expect_equal(WH_1d(d, ec, criterion = "AIC"),
             WH_1d_outer(d, ec, criterion = "AIC"))
expect_equal(WH_1d(d, ec, criterion = "BIC"),
             WH_1d_outer(d, ec, criterion = "BIC"))
expect_equal(WH_1d(d, ec, criterion = "GCV"),
             WH_1d_outer(d, ec, criterion = "GCV"))

## rank reduction works----
ref_perf_red <- WH_1d_perf(d, ec, p = 20)
ref_outer_red <- WH_1d_outer(d, ec, p = 20)
expect_equal(WH_1d(d, ec, p = 20), ref_outer_red)
expect_equal(WH_1d(d, ec, method = "perf", p = 20), ref_perf_red)
expect_equal(ref_perf_red, ref_perf_red, tolerance = 1e-1)

# Extrapolation----

newdata <- 18:99

## Regression----
perf_extra_reg <- WH_1d_perf(y = y, wt = wt, reg = TRUE) |> predict(newdata)
outer_extra_reg <- WH_1d_outer(y = y, wt = wt, reg = TRUE) |> predict(newdata)
expect_equal(perf_extra_reg, outer_extra_reg, tolerance = 1e-5)

## Maximum likelihood----
perf_extra_ml <- WH_1d_perf(d, ec) |> predict(newdata)
outer_extra_ml <- WH_1d_outer(d, ec) |> predict(newdata)
expect_equal(perf_extra_ml, outer_extra_ml, tolerance = 1e-1)

# Plots----

## Regression----
WH_1d_outer(y = y, wt = wt, reg = TRUE) |> plot()
WH_1d_perf(y = y, wt = wt, reg = TRUE) |> plot()

WH_1d_outer(y = y, wt = wt, reg = TRUE) |> plot("res")
WH_1d_perf(y = y, wt = wt, reg = TRUE) |> plot("res")

WH_1d_outer(y = y, wt = wt, reg = TRUE) |> plot("edf")
WH_1d_perf(y = y, wt = wt, reg = TRUE) |> plot("edf")

## Maximum likelihood----
WH_1d_outer(d, ec) |> plot()
WH_1d_perf(d, ec) |> plot()

WH_1d_outer(d, ec) |> plot("res")
WH_1d_perf(d, ec) |> plot("res")

WH_1d_outer(d, ec) |> plot("edf")
WH_1d_perf(d, ec) |> plot("edf")

## Extrapolation----
perf_extra_reg |> plot()
outer_extra_reg |> plot()

perf_extra_ml |> plot()
outer_extra_ml |> plot()

