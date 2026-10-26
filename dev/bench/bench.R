# Timings of WH on simulated 2D grids (Poisson deaths, log-hazard additive in both dimensions).
#   Rscript dev/bench/bench.R [library path]   (library path : where the WH under test is installed)
# Run each version on the same machine, in the same job, to compare : the machine is shared.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) .libPaths(c(args[[1]], .libPaths()))
library(WH)
sim2 <- function(n1, n2, seed = 1) {
  set.seed(seed)
  x <- seq_len(n1); z <- seq_len(n2)
  eta <- outer(-6 + 0.08 * x, 0.3 * exp(-z / 5), "+")
  ec <- matrix(runif(n1 * n2, 50, 500), n1, n2)
  d <- matrix(rpois(n1 * n2, ec * exp(eta)), n1, n2)
  dimnames(d) <- dimnames(ec) <- list(x, z); list(d = d, ec = ec)
}
el <- function(expr) system.time(expr)[["elapsed"]]
cat("WH", as.character(packageVersion("WH")), "\n")
for (sz in list(c(60, 30), c(100, 40))) {
  s <- sim2(sz[1], sz[2])
  t_fit <- el(o <- WH(s$d, s$ec, verbose = 0))
  t_dv <- el(for (i in 1:3) WH:::diag_V_compact_cpp(o$R)) / 3
  t_fix <- el(WH(s$d, s$ec, lambda = o$lambda, verbose = 0))
  t_pr <- el(op <- predict(o, newdata = list(seq_len(sz[1] + 10), seq_len(sz[2] + 5))))
  cat(sprintf("%dx%d p=%d band=%d evals=%d | outer %.2fs  fixed %.2fs  diag_V %.3fs  predict %.2fs\n",
      sz[1], sz[2], prod(sz), nrow(o$R) - 1, length(o$lambda_track), t_fit, t_fix, t_dv, t_pr))
}
