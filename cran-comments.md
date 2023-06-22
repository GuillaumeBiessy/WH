## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission following recent package acceptance but failure to additional tests using the alternative BLAS / LAPACK implementation from Intel MKL '2022.3.0'.

In the two previous versions of the package (1.0.3 and 1.0.4), failure of user-facing function tests were observed when using MKL. The failed tests are of the form "Test that calling user-facing function A which calls internal function B is the same as calling internal function B directly". As the functions of the `WH` package are fully deterministic, those tests should either fail because function A does not properly call function B, or pass as both objects make the same function call in the end. I made sure of that by temporarily replacing `test_that` by `test_identical` in the tests and those tests still passed.

However, when MKL is used, those tests fail and show a relative error close to 1e-7 on the results of those two similar calls. This is higher than the default tolerance of 1e-8 of function `test_that` and therefore it makes the test fail. Successive calls to the same function do not seem to return identical results when MKL is used. Hence, the only solution I see for making all tests still work for MKL is building tests that better deal with this variability in the results by having higher tolerance.

In this version I have:

* Further improved the robustness of the model tests by only testing the model fit, estimated standard deviation (with a 100 times higher relative tolerance) and REML criterion while comparing two fits. This should hopefully prevent test failure when using the MKL implementation of BLAS / LAPACK.

Please note that :

* Unfortunately I did not understand the comment that my tests lack a summary. Since version 1.0.4, those tests are embedded in test_that calls and start with a statement about what I am testing. I have looked for every occurrence of the substring "test" in "Writing R Extensions", reread all 3 chapters of "R Packages, 2nd edition" about Testing and looked at the tests implemented by the popular `ggplot2` package but was unable to find any clue about how to improve my test structure.

* I am still unable to figure out an easy way to install MKL on my Windows machine and link it with R.
