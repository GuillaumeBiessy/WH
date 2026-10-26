# Gap between the fits of the WH under test and the 2.0.0 reference fits (tests/testthat/fixtures).
#   Rscript dev/ref-2.0.0/compare_ref.R [library path]
# Prints, by fit and component, the largest relative gap max|a - b| / max|b| and n_eval of both.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) .libPaths(c(args[[1]], .libPaths()))
library(WH)
source(file.path("tests", "testthat", "helper-reference.R"))
ref <- readRDS(file.path("tests", "testthat", "fixtures", "ref-2.0.0.rds"))
cur <- make_reference_fits()
rel <- function(a, b) if (is.null(a) || is.null(b)) NA else max(abs(unlist(a) - unlist(b))) / max(abs(unlist(b)), 1e-300)
out <- do.call(rbind, lapply(names(ref), \(k) {
  comps <- c("lambda", "y_hat", "std_y_hat", "edf", "y_pred", "std_y_pred")
  gaps <- sapply(comps, \(cmp) rel(cur[[k]][[cmp]], ref[[k]][[cmp]]))
  crit <- rel(cur[[k]]$diagnosis[c("dev", "pen", "AIC", "BIC", "GCV", "REML")], ref[[k]]$diagnosis[c("dev", "pen", "AIC", "BIC", "GCV", "REML")])
  data.frame(fit = k, t(signif(gaps, 2)), diagnosis = signif(crit, 2), n_eval = paste(ref[[k]]$n_eval, cur[[k]]$n_eval, sep = "/"))
}))
cat("WH", as.character(packageVersion("WH")), "against the 2.0.0 reference\n")
print(out, row.names = FALSE)
