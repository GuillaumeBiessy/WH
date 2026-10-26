# Builds the fits compared with the WH 2.0.0 reference (fixtures/ref-2.0.0.rds): the same calls as
# dev/ref-2.0.0/make_ref.R, which saved the reference with the installed 2.0.0.
make_reference_fits <- function() {

  keep <- function(o) {
    out <- o[intersect(names(o), c("y_hat", "std_y_hat", "res", "edf", "lambda", "diagnosis", "y_pred", "std_y_pred"))]
    out$n_eval <- length(o$lambda_track)
    out
  }
  fit <- function(...) suppressWarnings(suppressMessages(WH(..., verbose = 0)))

  d1 <- portfolio_mort$d; ec1 <- portfolio_mort$ec
  y1 <- log(d1 / ec1); wt1 <- d1
  d2 <- portfolio_LTC$d; ec2 <- portfolio_LTC$ec
  y2 <- log(d2 / ec2); y2[d2 == 0] <- - 20; wt2 <- d2

  ref <- list()
  for (q in 1:3) for (crit in c("REML", "AIC", "BIC", "GCV")) {
    ref[[sprintf("1d_ml_q%d_%s", q, crit)]] <- keep(fit(d1, ec1, q = q, criterion = crit))
    ref[[sprintf("1d_reg_q%d_%s", q, crit)]] <- keep(fit(y = y1, wt = wt1, q = q, criterion = crit))
  }
  for (crit in c("REML", "GCV")) {
    ref[[sprintf("2d_ml_q22_%s", crit)]] <- keep(fit(d2, ec2, criterion = crit))
    ref[[sprintf("2d_reg_q22_%s", crit)]] <- keep(fit(y = y2, wt = wt2, criterion = crit))
  }
  ref[["2d_ml_q13_REML"]] <- keep(fit(d2, ec2, q = c(1, 3)))
  ref[["1d_ml_fixed"]] <- keep(fit(d1, ec1, lambda = 1e2))
  ref[["2d_ml_fixed"]] <- keep(fit(d2, ec2, lambda = c(1e1, 1e2)))
  ref[["1d_ml_pred"]] <- keep(predict(fit(d1, ec1), newdata = 40:99))
  ref[["2d_ml_pred"]] <- keep(predict(fit(d2, ec2), newdata = list(age = 60:109, duration = 0:19)))
  ref[["2d_reg_pred"]] <- keep(predict(fit(y = y2, wt = wt2), newdata = list(age = 60:109, duration = 0:19)))
  ref
}
