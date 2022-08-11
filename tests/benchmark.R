# data("portfolios_mort")
# data("portfolios_LTC")
#
# library(profvis)
# library(bench)
#
# # 1D smoothing----
# ec <- rowSums(portfolios_mort[[1]]$expo) / 365.25
# d <- rowSums(portfolios_mort[[1]]$exit)
#
# keep <- which(ec > 0) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()
# ec <- ec[keep]
# d <- d[keep]
#
# y <- log(d / ec)
# y[d == 0] <- - 20
# wt <- d
#
# mark(WH_1d(d, ec, framework = "reg"))
# mark(WH_1d(d, ec, framework = "reg", method = "optim"))
#
# mark(WH_1d(d, ec))
# mark(WH_1d(d, ec, method = "optim"))
#
# library(tidyverse)
#
# # Tuning of convergence parameters for regression
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_1d_reg_optim(
#     y = y, wt = wt, accu_edf = 10 ^ accu)),
#     lambda = map_dbl(fit, "lambda"),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_1d_reg_optim(
#       y = y, wt = wt, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit)
#
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_1d_reg_fs(
#     y = y, wt = wt, accu_edf = 10 ^ accu)),
#     lambda = map_dbl(fit, "lambda"),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_1d_reg_fs(
#       y = y, wt = wt, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit)
#
# # Tuning of convergence parameters for maximum likelihood
#
# ## edf
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_1d_ml_optim(
#     d = d, ec = ec, accu_edf = 10 ^ accu)),
#     lambda = map_dbl(fit, "lambda"),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_1d_ml_optim(
#       d = d, ec = ec, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit)
#
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_1d_ml_fs(
#     d = d, ec = ec, accu_edf = 10 ^ accu)),
#     lambda = map_dbl(fit, "lambda"),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_1d_ml_fs(
#       d = d, ec = ec, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit)
#
# ## dev
# tibble(accu_dev = (- 12):0) |>
#   mutate(fit = map(accu_dev, \(accu) WH_1d_ml_optim(
#     d = d, ec = ec, accu_dev = 10 ^ accu)),
#     lambda = map_dbl(fit, "lambda"),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_dev, \(accu) mark(WH_1d_ml_optim(
#       d = d, ec = ec, accu_dev = 10 ^ accu))$median[[1]])) |>
#   select(- fit)
#
# tibble(accu_dev = (- 12):0) |>
#   mutate(fit = map(accu_dev, \(accu) WH_1d_ml_fs(
#     d = d, ec = ec, accu_dev = 10 ^ accu)),
#     lambda = map_dbl(fit, "lambda"),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_dev, \(accu) mark(WH_1d_ml_fs(
#       d = d, ec = ec, accu_dev = 10 ^ accu))$median[[1]])) |>
#   select(- fit)
#
# # 2D smoothing----
# d  <- (portfolios_LTC[[1]]$exit) |> aperm(c(3,4,1,2)) |> colSums(dims = 2)
# ec <- (portfolios_LTC[[1]]$expo / 365.25) |> aperm(c(3,1,2)) |> colSums()
#
# keep_age <- which(rowSums(ec) > 1e2) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()
# keep_duration <- which(colSums(ec) > 1e2) |> range() |> (\(x) seq(x[[1]], x[[2]], 1))()
#
# d <- d[keep_age, keep_duration]
# ec <- ec[keep_age, keep_duration]
#
# y <- log(d / ec) # observation vector
# y[d == 0] <- - 20
# wt <- d
#
# mark(WH_2d(d, ec, framework = "reg"))
# mark(WH_2d(d, ec, framework = "reg", method = "optim"))
#
# mark(WH_2d(d, ec, framework = "reg", p = c(10, 5)))
# mark(WH_2d(d, ec, framework = "reg", method = "optim", p = c(10, 5)))
#
# mark(WH_2d(d, ec))
# mark(WH_2d(d, ec, method = "optim"))
#
# mark(WH_2d(d, ec, p = c(10, 5)))
# mark(WH_2d(d, ec, method = "optim", p = c(10, 5)))
#
# # Tuning of convergence parameters for regression
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_2d_reg_optim(
#     d = d, ec = ec, accu_edf = 10 ^ accu)),
#     lambda = map(fit, "lambda"),
#     lambda1 = map_dbl(lambda, 1),
#     lambda2 = map_dbl(lambda, 2),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_2d_reg_optim(
#       d = d, ec = ec, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit, - lambda)
#
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_2d_reg_fs(
#     d = d, ec = ec, accu_edf = 10 ^ accu)),
#     lambda = map(fit, "lambda"),
#     lambda1 = map_dbl(lambda, 1),
#     lambda2 = map_dbl(lambda, 2),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_2d_reg_fs(
#       d = d, ec = ec, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit, - lambda)
#
# # Tuning of convergence parameters for maximum likelihood
#
# ## edf
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_2d_ml_optim(
#     d = d, ec = ec, accu_edf = 10 ^ accu)),
#     lambda = map(fit, "lambda"),
#     lambda1 = map_dbl(lambda, 1),
#     lambda2 = map_dbl(lambda, 2),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_2d_ml_optim(
#       d = d, ec = ec, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit, - lambda)
#
# tibble(accu_edf = (- 12):0) |>
#   mutate(fit = map(accu_edf, \(accu) WH_2d_ml_fs(
#     d = d, ec = ec, accu_edf = 10 ^ accu)),
#     lambda = map(fit, "lambda"),
#     lambda1 = map_dbl(lambda, 1),
#     lambda2 = map_dbl(lambda, 2),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_edf, \(accu) mark(WH_2d_ml_fs(
#       d = d, ec = ec, accu_edf = 10 ^ accu))$median[[1]])) |>
#   select(- fit, - lambda)
#
# ## dev
# tibble(accu_dev = (- 12):0) |>
#   mutate(fit = map(accu_dev, \(accu) WH_2d_ml_optim(
#     d = d, ec = ec, accu_dev = 10 ^ accu)),
#     lambda = map(fit, "lambda"),
#     lambda1 = map_dbl(lambda, 1),
#     lambda2 = map_dbl(lambda, 2),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_dev, \(accu) mark(WH_2d_ml_optim(
#       d = d, ec = ec, accu_dev = 10 ^ accu))$median[[1]])) |>
#   select(- fit, - lambda)
#
# tibble(accu_dev = (- 12):0) |>
#   mutate(fit = map(accu_dev, \(accu) WH_2d_ml_fs(
#     d = d, ec = ec, accu_dev = 10 ^ accu)),
#     lambda = map(fit, "lambda"),
#     lambda1 = map_dbl(lambda, 1),
#     lambda2 = map_dbl(lambda, 2),
#     edf = map(fit, "edf") |> map_dbl(sum),
#     time_spent = map_dbl(accu_dev, \(accu) mark(WH_2d_ml_fs(
#       d = d, ec = ec, accu_dev = 10 ^ accu))$median[[1]])) |>
#   select(- fit, - lambda)

# # Plots----
# library(ggplot2)
# library(Linkplots)
#
# theme_LinkPact() |> theme_set()
#
# fit_names <- c("Régression / GCV", "Régression / REML", "Vraisemblance / GCV", "Vraisemblance / REML")
#
# # 1D
#
# list(WH_1d_fit_reg_optim_pred, WH_1d_fit_reg_fs_pred, WH_1d_fit_ml_optim_pred, WH_1d_fit_ml_fs_pred) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   tibble::tibble() |>
#   dplyr::mutate(y_hat_diff = (exp(y_pred_2) - exp(y_pred)) / exp(y_pred),
#                 std_y_hat_diff = (std_y_pred_2 - std_y_pred) / std_y_pred) |>
#   tidyr::pivot_longer(cols = c(y_hat_diff, std_y_hat_diff)) |>
#   ggplot(aes(x = x, y = value, color = name)) +
#   geom_line() +
#   facet_wrap(~fit)
#
# list(WH_1d_fit_reg_optim, WH_1d_fit_reg_fs, WH_1d_fit_ml_optim, WH_1d_fit_ml_fs) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_lines(x, res, fill = fit)
#
# list(WH_1d_fit_reg_optim, WH_1d_fit_reg_fs, WH_1d_fit_ml_optim, WH_1d_fit_ml_fs) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_lines(x, std_y_hat, fill = fit)
#
# # 2D
#
# list(WH_2d_fit_reg_optim_pred, WH_2d_fit_reg_fs_pred, WH_2d_fit_ml_optim_pred, WH_2d_fit_ml_fs_pred) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   tibble::tibble() |>
#   dplyr::filter(t == 0) |>
#   dplyr::mutate(y_hat_diff = (exp(y_hat_pred) - exp(y_hat)) / exp(y_hat),
#                 std_y_hat_diff = (std_y_hat_pred - std_y_hat) / std_y_hat) |>
#   tidyr::pivot_longer(cols = c(y_hat_diff, std_y_hat_diff)) |>
#   ggplot(aes(x = x, y = value, color = name)) +
#   geom_line() +
#   facet_wrap(~fit)
#
#
# list(WH_2d_fit_reg_optim, WH_2d_fit_reg_fs, WH_2d_fit_ml_optim, WH_2d_fit_ml_fs) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_2d_tile(t, x, res |> pmax(- 2.5) |> pmin(2.5), fit, divergent = TRUE, breaks_fill = seq(- 2.5, 2.5, 0.5))
#
# list(WH_2d_fit_reg_optim, WH_2d_fit_reg_fs, WH_2d_fit_ml_optim, WH_2d_fit_ml_fs) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_2d_contour(t, x, edf, fit)
#
# list(WH_2d_fit_reg_optim, WH_2d_fit_reg_fs, WH_2d_fit_ml_optim, WH_2d_fit_ml_fs) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_2d_contour(t, x, y_hat, fit)
#
# list(WH_2d_fit_reg_optim, WH_2d_fit_reg_fs, WH_2d_fit_ml_optim, WH_2d_fit_ml_fs) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_2d_contour(t, x, std_y_hat, fit)
#
# list(WH_2d_fit_reg_optim_pred, WH_2d_fit_reg_fs_pred, WH_2d_fit_ml_optim_pred, WH_2d_fit_ml_fs_pred) |>
#   map(output_to_df) |>
#   map2(fit_names, \(x,y) tibble::add_column(x, fit = y, .before = 1)) |>
#   do.call(what = dplyr::bind_rows) |>
#   dplyr::mutate(fit = factor(fit, levels = fit_names)) |>
#   lkp_2d_contour(t, x, y_hat, fit)
