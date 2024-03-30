## R CMD check results

0 errors | 0 warnings | 0 notes

This minor patch aims to fix an issue with a single failed test of the WH package, which occurs in the additional test MKL.

This failure is linked to a puzzling behaviour of the MKL BLAS which returns random results (i.e. not the same from one call to another) when performing eigendecomposition of matrices of size greater than 47.

A minimal reproducing example for this issue is :
D_mat <- diff(diag(47), differences = 2)
P  <- crossprod(D_mat)
testthat::expect_identical(eigen(P), eigen(P)) # works

D_mat <- diff(diag(48), differences = 2)
P  <- crossprod(D_mat)
testthat::expect_equal(eigen(P), eigen(P)) # fails

Package WH was already failing this particular MKL test when first accepted (and then reaccepted) to CRAN.

In this version, I have :

* Increased the tolerance of some of the tests from 1e-6 to 1e-5
