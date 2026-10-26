# Calibration of the non-regression tolerance: how much does WH 2.0.0 itself move when its inputs
# move by 1e-14 relative? Run from the package root against an installed WH 2.0.0:
#   Rscript dev/ref-2.0.0/calib_noise.R
library(WH); stopifnot(as.character(packageVersion("WH")) == "2.0.0")
portfolio_mort <- WH::portfolio_mort
portfolio_mort$ec <- portfolio_mort$ec * (1 + 1e-14)
source(file.path("tests", "testthat", "helper-reference.R"))
ref <- readRDS(file.path("tests", "testthat", "fixtures", "ref-2.0.0.rds"))
cur <- make_reference_fits()
rel <- function(a, b) max(abs(unlist(a) - unlist(b))) / max(abs(unlist(b)))
for (k in grep("^1d", names(ref), value = TRUE)) cat(sprintf("%-16s lambda %.1e  y_hat %.1e  std_y_hat %.1e\n", k,
  rel(cur[[k]]$lambda, ref[[k]]$lambda), rel(cur[[k]]$y_hat, ref[[k]]$y_hat), rel(cur[[k]]$std_y_hat, ref[[k]]$std_y_hat)))
