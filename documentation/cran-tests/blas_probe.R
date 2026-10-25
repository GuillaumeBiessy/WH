# blas_probe.R -- prints, in one process started like the check's subprocesses, which BLAS / LAPACK answer.
# Run with MKL_VERBOSE=1 or OPENBLAS_VERBOSE=2 set: the library then prints its own banner on stdout.
x <- crossprod(matrix(seq_len(400) / 7, 20))
invisible(chol(x))
invisible(eigen(x, symmetric = TRUE))
cat("PROBE R_HOME=", R.home(), " LAPACK=", La_version(), "\n", sep = "")
