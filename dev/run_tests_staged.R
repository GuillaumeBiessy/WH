# Gate of the WH test suite, detected by the landing launcher (gate-lock-r, land.sh).
#
# Loads the package from this tree with pkgload::load_all, runs every
# tests/testthat/test-*.R and counts failures AND errors.
#
# Exit contract (fail-closed): status 0 only if the verdict starts with "green".
# "green" requires zero failure, zero error and zero skipped test: a green run that
# skips tests does not prove the coverage it claims. The verdict is also written
# atomically to dev/.last_test_status, durations by file to dev/.test_durations.tsv.
#
# Usage, from the package root (under the gate-lock-r lock when run in background):
#   Rscript dev/run_tests_staged.R > dev/.last_test_run.log 2>&1

if (!file.exists("DESCRIPTION") || read.dcf("DESCRIPTION")[1, "Package"] != "WH") {
  stop("Run from the root of the WH package.")
}
root <- normalizePath(".", winslash = "/")
status_file <- file.path(root, "dev", ".last_test_status")
durations_file <- file.path(root, "dev", ".test_durations.tsv")

write_status <- function(line) {
  tmp <- paste0(status_file, ".tmp")
  writeLines(line, tmp)
  invisible(file.rename(tmp, status_file))
}
write_status(paste("STARTED @", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

verdict <- tryCatch({
  # Optimised build (load_all alone compiles in debug mode, -O0: the 2D tests ran ten times slower)
  pkgbuild::compile_dll(root, debug = FALSE, quiet = TRUE)
  suppressMessages(pkgload::load_all(root, export_all = TRUE, helpers = FALSE, quiet = TRUE, compile = FALSE))
  loaded_path <- normalizePath(getNamespaceInfo(asNamespace("WH"), "path"), winslash = "/")
  if (!identical(tolower(loaded_path), tolower(root))) {
    stop("Package loaded from ", loaded_path, " instead of ", root)
  }

  files <- list.files(file.path(root, "tests", "testthat"), pattern = "^test-.*\\.R$", full.names = TRUE)
  if (length(files) == 0) stop("No test file under tests/testthat/.")
  totals <- c(pass = 0, fail = 0, skip = 0)
  durations <- data.frame(file = character(), seconds = numeric())
  for (f in files) {
    t0 <- Sys.time()
    res <- tryCatch(
      as.data.frame(testthat::test_file(f, reporter = "summary", stop_on_failure = FALSE)),
      error = function(err) {
        message("ERROR outside a test block in ", basename(f), ": ", conditionMessage(err))
        data.frame(nb = 1, failed = 1, error = TRUE, skipped = FALSE)
      }
    )
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    n_fail <- sum(res$failed) + sum(as.logical(res$error))
    n_skip <- sum(as.logical(res$skipped))
    totals <- totals + c(sum(res$nb) - n_fail - n_skip, n_fail, n_skip)
    durations[nrow(durations) + 1, ] <- list(basename(f), round(secs, 1))
    cat(sprintf("FILE %-28s fail=%d skip=%d secs=%.1f\n", basename(f), n_fail, n_skip, secs))
  }
  write.table(durations, durations_file, sep = "\t", row.names = FALSE, quote = FALSE)

  sprintf("%s pass=%d fail=%d skip=%d",
          if (totals[["fail"]] == 0 && totals[["skip"]] == 0) "green" else "red",
          totals[["pass"]], totals[["fail"]], totals[["skip"]])
}, error = function(err) paste("red error:", conditionMessage(err)))

write_status(verdict)
cat(verdict, "\n")
if (!interactive()) quit(status = if (startsWith(verdict, "green ")) 0 else 1)
