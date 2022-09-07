data("portfolio_LTC")
if(!interactive()) pdf(NULL)

# Data----
keep_age <- which(rowSums(portfolio_LTC$ec) > 1e2)
keep_duration <- which(colSums(portfolio_LTC$ec) > 1e2)

d  <- portfolio_LTC$d[keep_age, keep_duration]
ec <- portfolio_LTC$ec[keep_age, keep_duration]

y <- log(d / ec) # observation vector
y[d == 0] <- - 20
wt <- d

# Regression----

## Various way of calling regresion work and method with fixed lambda as well----
ref_fixed_lambda <- WH_2d_fixed_lambda(y = y, wt = wt, lambda = c(1e2,1e2), reg = TRUE)
expect_equal(WH_2d(y = y, wt = wt, lambda = c(1e2,1e2)), ref_fixed_lambda)
expect_equal(WH_2d(d, ec, framework = "reg", lambda = c(1e2,1e2)), ref_fixed_lambda)
expect_equal(WH_2d(d, y = y, lambda = c(1e2,1e2)), ref_fixed_lambda)
expect_equal(WH_2d_fixed_lambda(y = y, wt = wt, lambda = c(1e2,1e2), reg = TRUE), ref_fixed_lambda)

## FS is the default method and call FS----
ref_fs <- WH_2d_fs(y = y, wt = wt, reg = TRUE)
expect_equal(WH_2d(y = y, wt = wt, method = "fs"), ref_fs)
expect_equal(WH_2d(y = y, wt = wt), ref_fs)

## optim method call optim----
ref_outer <- WH_2d_outer(y = y, wt = wt, reg = TRUE)
expect_equal(WH_2d(y = y, wt = wt, method = "outer"), ref_outer)

## optim and fs method match for regression----
expect_equal(ref_fs, ref_outer, tolerance = 1e-4)

## REML is default criterion for optim----
expect_equal(WH_2d(y = y, wt = wt, method = "outer", criterion = "REML"), ref_outer)

## other criteria work----
expect_equal(WH_2d(y = y, wt = wt, criterion = "AIC"),
             WH_2d_perf(y = y, wt = wt, criterion = "AIC", reg = TRUE))
expect_equal(WH_2d(y = y, wt = wt, criterion = "BIC"),
             WH_2d_perf(y = y, wt = wt, criterion = "BIC", reg = TRUE))
expect_equal(WH_2d(y = y, wt = wt, criterion = "GCV"),
             WH_2d_perf(y = y, wt = wt, criterion = "GCV", reg = TRUE))

## rank reduction works----
ref_fs_red <- WH_2d_fs(y = y, wt = wt, p = c(10, 5), reg = TRUE)
ref_outer_red <- WH_2d_outer(y = y, wt = wt, p = c(10, 5), reg = TRUE)

expect_equal(WH_2d(y = y, wt = wt, p = c(10, 5)), ref_fs_red)
expect_equal(WH_2d(y = y, wt = wt, method = "outer", p = c(10, 5)), ref_outer_red)
expect_equal(ref_fs_red, ref_outer_red, tolerance = 1e-4)

expect_equal(WH_2d(y = y, wt = wt, method = "fs", max_dim = 100),
             WH_2d(y = y, wt = wt, method = "outer", max_dim = 100), tolerance = 1e-4)

# Maximum likelihood----

## fixed lambda method works----
expect_equal(WH_2d(d, ec, lambda = c(1e2,1e2)),
             WH_2d_fixed_lambda(d, ec, lambda = c(1e2,1e2)))

## FS method with rank reduction works----
ref_fs_red <- WH_2d_fs(d, ec, p = c(10, 5))
expect_equal(WH_2d(d, ec, p = c(10, 5)), ref_fs_red)

## optim method with rank reduction works----
ref_outer_red <- WH_2d_outer(d, ec, p = c(10, 5))
expect_equal(WH_2d(d, ec, method = "outer", p = c(10, 5)), ref_outer_red)

## optim and fs method are not too far away for ML----
expect_equal(ref_fs_red, ref_outer_red, tolerance = 2e-1)

## automatic rank reduction works----
rr_fs <- WH_2d(d, ec, max_dim = 100)
rr_outer <- WH_2d(d, ec, method = "outer", max_dim = 100)
expect_equal(rr_fs, rr_outer, tolerance = 2e-1)

# Extrapolation----

newdata <- list(age = 50:99, duration = 0:19)

extra_fs <- rr_fs |> predict(newdata)
extra_outer <- rr_outer |> predict(newdata)

expect_equal(extra_fs,
             extra_outer,
             tolerance = 2e-1)

expect_equal(rr_fs |> predict(newdata),
             rr_fs |> predict_WH_2d_old(newdata), tolerance = 1e-2)

expect_equal(rr_fs |> predict(newdata, n_coef = 10),
             rr_fs |> predict_WH_2d_old(newdata), tolerance = 5e-3)
# Plots----

rr_fs |> plot()
rr_fs |> plot("std_y_hat")

extra_fs |> plot()
extra_fs |> plot("std_y_hat")

ref_fs |> predict(newdata) |> plot()
ref_fs |> predict(newdata) |> plot()


