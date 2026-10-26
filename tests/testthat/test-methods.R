# S3 methods and output formatting.

if (!interactive()) pdf(NULL)

fit1 <- WH(portfolio_mort$d, portfolio_mort$ec, verbose = 0)
fit2 <- WH(portfolio_LTC$d, portfolio_LTC$ec, verbose = 0)
pred1 <- predict(fit1, newdata = 40:99)
pred2 <- predict(fit2, newdata = list(age = 60:109, duration = 0:19))

test_that("print summarises the fit", {

  expect_output(print(fit1), "Initial data contains 45 data points")
  expect_output(print(fit1), "Observation positions:  50  to  94")
  expect_output(print(fit2), "Initial data contains 450 data points")
  expect_invisible(print(fit2))
})

test_that("every plot type draws without error, before and after prediction", {

  for (what in c("fit", "res", "edf")) expect_no_error(plot(fit1, what))
  expect_no_error(plot(pred1))
  expect_no_error(plot(fit1, trans = exp))
  for (what in c("y_hat", "std_y_hat", "res", "edf")) expect_no_error(plot(fit2, what))
  expect_no_error(plot(pred2))
  expect_error(plot(fit1, "y_hat"), "what should be one of")
})

test_that("output_to_df gathers fitted and predicted values", {

  df1 <- output_to_df(fit1)
  expect_equal(nrow(df1), 45)
  expect_equal(df1$y_hat, unname(c(fit1$y_hat)))
  dfp <- output_to_df(pred1)
  expect_equal(nrow(dfp), 60)
  expect_true(all(is.na(dfp$y[dfp$x < 50])))
  df2 <- output_to_df(pred2, dim1 = "age", dim2 = "duration")
  expect_equal(nrow(df2), 50 * 20)
  expect_named(df2, c("age", "duration", "y", "y_hat", "std_y_hat", "res", "edf"))
})

test_that("vcov matches the standard deviations, with and without prediction", {

  expect_equal(sqrt(diag(vcov(fit1))), unname(fit1$std_y_hat), tolerance = 1e-8)
  expect_equal(sqrt(diag(vcov(pred1))), unname(pred1$std_y_pred), tolerance = 1e-8)
  expect_equal(vcov(pred1, pred = FALSE), vcov(fit1))
  expect_equal(sqrt(diag(vcov(pred2))), c(pred2$std_y_pred), tolerance = 1e-8)
  expect_equal(vcov(pred2, pred = FALSE), vcov(fit2))
})
