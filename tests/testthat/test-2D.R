data("portfolios_LTC")

# 2D smoothing----
keep_age <- which(rowSums(portfolios_LTC[[1]]$ec) > 1e2)
keep_duration <- which(colSums(portfolios_LTC[[1]]$ec) > 1e2)

d  <- portfolios_LTC[[1]]$d[keep_age, keep_duration]
ec <- portfolios_LTC[[1]]$ec[keep_age, keep_duration]

y <- log(d / ec) # observation vector
y[d == 0] <- - 20
wt <- d

# Regression

ref_fixed_lambda <- WH_2d_reg_fixed_lambda(y, wt, lambda = c(1e2,1e2))
expect_equal(WH_2d(y = y, wt = wt, lambda = c(1e2,1e2)), ref_fixed_lambda)
expect_equal(WH_2d(d, ec, framework = "reg", lambda = c(1e2,1e2)), ref_fixed_lambda)
expect_equal(WH_2d(d, y = y, lambda = c(1e2,1e2)), ref_fixed_lambda)

ref_fs <- WH_2d_reg_fs(y, wt)
expect_equal(WH_2d(y = y, wt = wt, method = "fs"), ref_fs)
expect_equal(WH_2d(y = y, wt = wt), ref_fs)

ref_optim <- WH_2d_reg_optim(y, wt)
expect_equal(WH_2d(y = y, wt = wt, method = "optim"), ref_optim)

expect_equal(ref_fs, ref_optim, tolerance = 1e-4)

expect_equal(WH_2d(y = y, wt = wt, method = "optim", criterion = "REML"), ref_optim)

ref_aic <- WH_2d_reg_optim(y, wt, criterion = "AIC")
expect_equal(WH_2d(y = y, wt = wt, criterion = "AIC"), ref_aic)

ref_bic <- WH_2d_reg_optim(y, wt, criterion = "BIC")
expect_equal(WH_2d(y = y, wt = wt, criterion = "BIC"), ref_bic)

ref_gcv <- WH_2d_reg_optim(y, wt, criterion = "GCV")
expect_equal(WH_2d(y = y, wt = wt, criterion = "GCV"), ref_gcv)

ref_fs_red <- WH_2d_reg_fs(y, wt, p = c(10, 5))
expect_equal(WH_2d(y = y, wt = wt, p = c(10, 5)), ref_fs_red)

ref_optim_red <- WH_2d_reg_optim(y, wt, p = c(10, 5))
expect_equal(WH_2d(y = y, wt = wt, method = "optim", p = c(10, 5)), ref_optim_red)

expect_equal(ref_fs_red, ref_optim_red, tolerance = 1e-4)

# Maximum likelihood
expect_equal(WH_2d(d, ec, lambda = c(1e2,1e2)),
             WH_2d_ml_fixed_lambda(d, ec, lambda = c(1e2,1e2)))

ref_fs_red <- WH_2d_ml_fs(d, ec, p = c(10, 5))
expect_equal(WH_2d(d, ec, p = c(10, 5)), ref_fs_red)

ref_optim_red <- WH_2d_ml_optim(d, ec, p = c(10, 5))
expect_equal(WH_2d(d, ec, method = "optim", p = c(10, 5)), ref_optim_red)

expect_equal(ref_fs_red, ref_optim_red, tolerance = 1e-1)

# Plots
WH_2d(d, ec) |> plot()
WH_2d(d, ec) |> plot("std_y_hat")

# Extrapolation
newdata <- list(age = 50:99, duration = 0:19)

extra_fs <- WH_2d(d, ec) |> predict(newdata)
extra_optim <- WH_2d(d, ec, method = "optim") |> predict(newdata)

expect_equal(extra_fs,
             extra_optim,
             tolerance = 2e-1)

plot(extra_fs)
plot(extra_optim)

# Rank reduction
keep_age <- which(rowSums(portfolios_LTC[[1]]$ec) > 0)
keep_duration <- which(colSums(portfolios_LTC[[1]]$ec) > 0)

d  <- portfolios_LTC[[1]]$d[keep_age, keep_duration]
ec <- portfolios_LTC[[1]]$ec[keep_age, keep_duration]

prod(dim(d)) # dimension of problem is 1232

expect_equal(WH_2d(y = y, wt = wt, method = "fs", max_dim = 100),
             WH_2d(y = y, wt = wt, method = "optim", max_dim = 100), tolerance = 1e-4)
expect_equal(WH_2d(d, ec, max_dim = 100), WH_2d(d, ec, method = "optim", max_dim = 100), tolerance = 2e-1)
