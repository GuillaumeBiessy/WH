# Gaps on the total effective degrees of freedom sum(edf) against the 2.0.0 reference fits
# (tests/testthat/fixtures), the yardstick of the selected-lambda tolerance of test-reference.R.
#   Rscript dev/ref-2.0.0/edf_gaps.R [library path]   WH under test, inputs as they are
#   Rscript dev/ref-2.0.0/edf_gaps.R --calib          installed WH 2.0.0, inputs perturbed by 1e-14
# --calib also runs 2.0.0 on unperturbed inputs, which must reproduce the reference exactly.
# Perturbations: exposures scaled by 1 +/- 1e-14 (a shift of the offset, which leaves the fit
# unchanged in exact arithmetic) and exposures times 1 + 1e-14 U(-1, 1) (two seeds), on both
# portfolios. The largest gap over the perturbations is reported.
args <- commandArgs(trailingOnly = TRUE)
calib <- identical(args, "--calib")
if (length(args) == 1 && !calib) .libPaths(c(args[[1]], .libPaths()))
library(WH)
if (calib) stopifnot(as.character(packageVersion("WH")) == "2.0.0")
source(file.path("tests", "testthat", "helper-reference.R"))
ref <- readRDS(file.path("tests", "testthat", "fixtures", "ref-2.0.0.rds"))

gap_edf <- function(cur) sapply(names(ref), \(k) abs(sum(cur[[k]]$edf) - sum(ref[[k]]$edf)))
gap_rel <- function(cur, cmp) sapply(names(ref), \(k) max(abs(unlist(cur[[k]][[cmp]]) - unlist(ref[[k]][[cmp]]))) / max(abs(unlist(ref[[k]][[cmp]]))))

with_inputs <- function(f_ec) {
  for (nm in c("portfolio_mort", "portfolio_LTC")) {
    data <- getExportedValue("WH", nm); data$ec <- f_ec(data$ec)
    assign(nm, data, envir = globalenv())
  }
  make_reference_fits()
}

if (calib) {
  exact <- with_inputs(identity)
  stopifnot(max(gap_edf(exact)) == 0)
  perturb <- list(
    up = \(ec) ec * (1 + 1e-14),
    down = \(ec) ec * (1 - 1e-14),
    noise1 = \(ec) { set.seed(1); ec * (1 + 1e-14 * runif(length(ec), -1, 1)) },
    noise2 = \(ec) { set.seed(2); ec * (1 + 1e-14 * runif(length(ec), -1, 1)) })
  runs <- lapply(perturb, with_inputs)
  edf_gap <- do.call(pmax, lapply(runs, gap_edf))
  lambda_gap <- do.call(pmax, lapply(runs, gap_rel, "lambda"))
  cat("WH 2.0.0, inputs perturbed by 1e-14 (max over", length(runs), "perturbations), against the reference\n")
} else {
  cur <- with_inputs(identity)
  edf_gap <- gap_edf(cur)
  lambda_gap <- gap_rel(cur, "lambda")
  cat("WH", as.character(packageVersion("WH")), "against the 2.0.0 reference\n")
}
out <- data.frame(fit = names(ref), edf_total = signif(sapply(ref, \(r) sum(r$edf)), 4),
  gap_edf_total = signif(edf_gap, 2), gap_lambda_rel = signif(lambda_gap, 2))
print(out, row.names = FALSE)
cat(sprintf("largest gap on sum(edf): %.2e degree of freedom (%s)\n", max(edf_gap), names(ref)[which.max(edf_gap)]))
