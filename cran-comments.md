## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following recent package acceptance but failure to additional tests using the alternative BLAS / LAPACK implementation from Intel MKL '2022.3.0'.

In the two previous versions of the package (1.0.3 and 1.0.4), test failure was observed when using MKL. Some failed tests are of the form "Test that calling user-facing function A which calls internal function B is the same as calling internal function B directly". As the functions of the `WH` package are fully deterministic, those tests should either fail because function A does not properly call function B, or pass as both objects ultimately make the same call to function B in the end. I made sure of that by temporarily replacing `test_that` by `test_identical` in those tests which still passed.

However, when MKL is used, those tests nows fails because the relative error between the results of those two calls is higher than the default tolerance of 1e-8 of function `test_that`. Therefore it seems that successive calls to the same function may no longer be expected to return identical results when MKL is used. Hence, for making all tests work for MKL, I need to further improve the tolerance of the tests.

In this version I have:

* Further improved the robustness of the model tests by restricting testing to the fitted model value, the estimated standard deviation (with a 100 times higher relative tolerance) and the REML criterion while comparing two fits. This should hopefully prevent test failure when using the MKL implementation of BLAS / LAPACK.
